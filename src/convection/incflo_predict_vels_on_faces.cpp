#include <incflo.H>
#include <incflo_convection_K.H>

using namespace amrex;

void incflo::predict_vels_on_faces (int lev, Real time, MultiFab& u_mac, MultiFab& v_mac,
                                    MultiFab& w_mac, MultiFab& vel)
{
    // fillpatch on vel.
    const int ng = 2;
    fillpatch_velocity(lev, time, vel, ng);

#ifdef AMREX_USE_EB
    auto const& fact = this->EBFactory(lev);
    auto const& flags = fact.getMultiEBCellFlagFab();
    auto const& fcent = fact.getFaceCent();
#endif

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    {
        for (MFIter mfi(vel, TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            Box const& bx = mfi.tilebox();
            Box const& ubx = mfi.nodaltilebox(0);
            Box const& vbx = mfi.nodaltilebox(1);
            Box const& wbx = mfi.nodaltilebox(2);
            Array4<Real> const& u = u_mac.array(mfi);
            Array4<Real> const& v = v_mac.array(mfi);
            Array4<Real> const& w = w_mac.array(mfi);
            Array4<Real const> const& vcc = vel.const_array(mfi);
#ifdef AMREX_USE_EB
            EBCellFlagFab const& flagfab = flags[mfi];
            Array4<EBCellFlag const> const& flagarr = flagfab.const_array();
            if (flagfab.getType(bx) == FabType::covered)
            {
                amrex::ParallelFor(ubx, vbx, wbx,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept { u(i,j,k) = 0.0; },
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept { v(i,j,k) = 0.0; },
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept { w(i,j,k) = 0.0; });
            }
            else if (flagfab.getType(amrex::grow(bx,1)) == FabType::singlevalued)
            {
                Array4<Real const> const& fcx = fcent[0]->const_array(mfi);
                Array4<Real const> const& fcy = fcent[1]->const_array(mfi);
                Array4<Real const> const& fcz = fcent[2]->const_array(mfi);
                ::incflo_predict_vels_on_faces_eb(bx,ubx,vbx,wbx,u,v,w,vcc,
                                                  flagarr,fcx,fcy,fcz);
            }
            else
#endif
            {
                ::incflo_predict_vels_on_faces(ubx,vbx,wbx,u,v,w,vcc);
            }
        }
    }

    EB_set_covered_faces({&u_mac,&v_mac,&w_mac}, 0.0);
    VisMF::Write(u_mac, "u_mac");
    VisMF::Write(v_mac, "v_mac");
    VisMF::Write(w_mac, "w_mac");
    amrex::Abort("xxxxx predict_vels_on_faces: end of fillpatch_velocity");
}

void
incflo::incflo_predict_vels_on_faces ( int lev, Real time,
                                       Vector< std::unique_ptr<MultiFab> >& vel_in)
{
    BL_PROFILE("incflo::incflo_predict_vels_on_faces");

    Box domain(geom[lev].Domain());   

    iMultiFab cc_mask(grids[lev], dmap[lev], 1, 1);

    const int covered_value = 1;
    const int notcovered_value = 0;
    const int physical_boundaries_value = 0;
    const int interior_value = 1;

    cc_mask.BuildMask(geom[lev].Domain(), geom[lev].periodicity(),
        covered_value, notcovered_value,
        physical_boundaries_value, interior_value);

#ifdef AMREX_USE_EB
       // Get EB geometric info
       Array< const MultiCutFab*,AMREX_SPACEDIM> areafrac;
       Array< const MultiCutFab*,AMREX_SPACEDIM> facecent;

       areafrac  =   ebfactory[lev] -> getAreaFrac();
       facecent  =   ebfactory[lev] -> getFaceCent();
#endif

       Real small_vel = 1.e-10;
       Real  huge_vel = 1.2345e300;

       // ****************************************************************************
       // We will store the left and right states in arrays for interpolation
       // ****************************************************************************

#ifdef AMREX_USE_EB
       MultiFab upls(m_u_mac[lev]->boxArray(), dmap[lev], 1, 1, MFInfo(), *ebfactory[lev]);
       MultiFab umns(m_u_mac[lev]->boxArray(), dmap[lev], 1, 1, MFInfo(), *ebfactory[lev]);

       MultiFab vpls(m_v_mac[lev]->boxArray(), dmap[lev], 1, 1, MFInfo(), *ebfactory[lev]);
       MultiFab vmns(m_v_mac[lev]->boxArray(), dmap[lev], 1, 1, MFInfo(), *ebfactory[lev]);

       MultiFab wpls(m_w_mac[lev]->boxArray(), dmap[lev], 1, 1, MFInfo(), *ebfactory[lev]);
       MultiFab wmns(m_w_mac[lev]->boxArray(), dmap[lev], 1, 1, MFInfo(), *ebfactory[lev]);
#else
       MultiFab upls(m_u_mac[lev]->boxArray(), dmap[lev], 1, 1, MFInfo());
       MultiFab umns(m_u_mac[lev]->boxArray(), dmap[lev], 1, 1, MFInfo());

       MultiFab vpls(m_v_mac[lev]->boxArray(), dmap[lev], 1, 1, MFInfo());
       MultiFab vmns(m_v_mac[lev]->boxArray(), dmap[lev], 1, 1, MFInfo());

       MultiFab wpls(m_w_mac[lev]->boxArray(), dmap[lev], 1, 1, MFInfo());
       MultiFab wmns(m_w_mac[lev]->boxArray(), dmap[lev], 1, 1, MFInfo());
#endif

       // We need this just to avoid FPE (eg for outflow faces)
       upls.setVal(covered_val);
       umns.setVal(covered_val);
       vpls.setVal(covered_val);
       vmns.setVal(covered_val);
       wpls.setVal(covered_val);
       wmns.setVal(covered_val);

       // ****************************************************************************
       // First compute the slopes
       // ****************************************************************************
       int slopes_comp = 0; int num_comp = 3;
       incflo_compute_slopes(lev, time, *vel_in[lev], xslopes_u, yslopes_u, zslopes_u, 
                             slopes_comp, num_comp);

       // ****************************************************************************
       // Then predict to face centers
       // ****************************************************************************
       for (MFIter mfi(*vel_in[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
       {
          // Tilebox
          Box  bx = mfi.tilebox();

          Box ubx = mfi.nodaltilebox(0);
          Box vbx = mfi.nodaltilebox(1);
          Box wbx = mfi.nodaltilebox(2);

          Box ubx_grown = mfi.growntilebox(IntVect::TheDimensionVector(0));
          Box vbx_grown = mfi.growntilebox(IntVect::TheDimensionVector(1));
          Box wbx_grown = mfi.growntilebox(IntVect::TheDimensionVector(2));

#ifdef AMREX_USE_EB
          const EBFArrayBox&  vel_fab = static_cast<EBFArrayBox const&>((*vel_in[lev])[mfi]);
          const EBCellFlagFab&  flags = vel_fab.getEBCellFlagFab();

          if (flags.getType(amrex::grow(bx,0)) == FabType::covered )
          {
             Real val = 1.2345e300;
             const auto& umac_array = (*m_u_mac[lev])[mfi].array();
             const auto& vmac_array = (*m_v_mac[lev])[mfi].array();
             const auto& wmac_array = (*m_w_mac[lev])[mfi].array();
 
             AMREX_FOR_3D(ubx, i, j, k, { umac_array(i,j,k) = val; });
             AMREX_FOR_3D(vbx, i, j, k, { vmac_array(i,j,k) = val; });
             AMREX_FOR_3D(wbx, i, j, k, { wmac_array(i,j,k) = val; });
 
//           (*m_u_mac[lev])[mfi].setVal( 1.2345e300, ubx, 0, 1);
//           (*m_v_mac[lev])[mfi].setVal( 1.2345e300, vbx, 0, 1);
//           (*m_w_mac[lev])[mfi].setVal( 1.2345e300, wbx, 0, 1);
          }
  
          // No cut cells in this FAB
          else if (flags.getType(amrex::grow(bx,1)) == FabType::regular )
          {
#endif
             // Cell-centered velocity
             const auto& ccvel_fab = vel_in[lev]->array(mfi);

             // Cell-centered slopes
             const auto& xslopes_fab = (xslopes_u[lev])->array(mfi);
             const auto& yslopes_fab = (yslopes_u[lev])->array(mfi);
             const auto& zslopes_fab = (zslopes_u[lev])->array(mfi);

             // Face-centered velocity components
             const auto& umac_fab = (m_u_mac[lev])->array(mfi);
             const auto& vmac_fab = (m_v_mac[lev])->array(mfi);
             const auto& wmac_fab = (m_w_mac[lev])->array(mfi);

             // Face-centered left and right states
             const auto& upls_fab = upls.array(mfi);
             const auto& vpls_fab = vpls.array(mfi);
             const auto& wpls_fab = wpls.array(mfi);
             const auto& umns_fab = umns.array(mfi);
             const auto& vmns_fab = vmns.array(mfi);
             const auto& wmns_fab = wmns.array(mfi);

             // No cut cells in tile + 1-cell witdh halo -> use non-eb routine

             amrex::ParallelFor(ubx,[small_vel,ccvel_fab,xslopes_fab,upls_fab,umns_fab,umac_fab] 
               AMREX_GPU_DEVICE (int i, int j, int k) noexcept
             {
                 // X-faces
                 upls_fab(i,j,k) = ccvel_fab(i  ,j,k,0) - 0.5 * xslopes_fab(i  ,j,k,0);
                 umns_fab(i,j,k) = ccvel_fab(i-1,j,k,0) + 0.5 * xslopes_fab(i-1,j,k,0);
                 if ( umns_fab(i,j,k) < 0.0 && upls_fab(i,j,k) > 0.0 ) {
                    umac_fab(i,j,k) = 0.0;
                 } else {
                    Real avg = 0.5 * ( upls_fab(i,j,k) + umns_fab(i,j,k) );
                    if ( std::abs(avg) <  small_vel) { umac_fab(i,j,k) = 0.0;
                    } else if (avg >= 0)             { umac_fab(i,j,k) = umns_fab(i,j,k);
                    } else                           { umac_fab(i,j,k) = upls_fab(i,j,k);
                    }
                 }
             });

             amrex::ParallelFor(vbx,[small_vel,ccvel_fab,yslopes_fab,vpls_fab,vmns_fab,vmac_fab] 
               AMREX_GPU_DEVICE (int i, int j, int k) noexcept
             {
                 // Y-faces
                 vpls_fab(i,j,k) = ccvel_fab(i,j  ,k,1) - 0.5 * yslopes_fab(i,j  ,k,1);
                 vmns_fab(i,j,k) = ccvel_fab(i,j-1,k,1) + 0.5 * yslopes_fab(i,j-1,k,1);
                 if ( vmns_fab(i,j,k) < 0.0 && vpls_fab(i,j,k) > 0.0 ) {
                    vmac_fab(i,j,k) = 0.0;
                 } else {
                    Real avg = 0.5 * ( vpls_fab(i,j,k) + vmns_fab(i,j,k) );
                    if ( std::abs(avg) <  small_vel) { vmac_fab(i,j,k) = 0.0;
                    } else if (avg >= 0)             { vmac_fab(i,j,k) = vmns_fab(i,j,k);
                    } else                           { vmac_fab(i,j,k) = vpls_fab(i,j,k);
                    }
                 }
             });

             amrex::ParallelFor(wbx,[small_vel,ccvel_fab,zslopes_fab,wpls_fab,wmns_fab,wmac_fab] 
               AMREX_GPU_DEVICE (int i, int j, int k) noexcept
             {
                 // Z-faces
                 wpls_fab(i,j,k) = ccvel_fab(i,j,k  ,2) - 0.5 * zslopes_fab(i,j,k  ,2);
                 wmns_fab(i,j,k) = ccvel_fab(i,j,k-1,2) + 0.5 * zslopes_fab(i,j,k-1,2);
                 if ( wmns_fab(i,j,k) < 0.0 && wpls_fab(i,j,k) > 0.0 ) {
                    wmac_fab(i,j,k) = 0.0;
                 } else {
                    Real avg = 0.5 * ( wpls_fab(i,j,k) + wmns_fab(i,j,k) );
                    if ( std::abs(avg) <  small_vel) { wmac_fab(i,j,k) = 0.0;
                    } else if (avg >= 0)             { wmac_fab(i,j,k) = wmns_fab(i,j,k);
                    } else                           { wmac_fab(i,j,k) = wpls_fab(i,j,k);
                    }
                 }
             });
             
#ifdef AMREX_USE_EB
          // Cut cells in this FAB
          } else {

             const auto& ccvel_fab = vel_in[lev]->array(mfi);

             // Cell-centered slopes
             const auto& xslopes_fab = (xslopes_u[lev])->array(mfi);
             const auto& yslopes_fab = (yslopes_u[lev])->array(mfi);
             const auto& zslopes_fab = (zslopes_u[lev])->array(mfi);

             // Face-centered left and right states
             const auto& upls_fab = upls.array(mfi);
             const auto& vpls_fab = vpls.array(mfi);
             const auto& wpls_fab = wpls.array(mfi);
             const auto& umns_fab = umns.array(mfi);
             const auto& vmns_fab = vmns.array(mfi);
             const auto& wmns_fab = wmns.array(mfi);

             // Face-centered areas
             const auto& apx_fab = areafrac[0]->array(mfi);
             const auto& apy_fab = areafrac[1]->array(mfi);
             const auto& apz_fab = areafrac[2]->array(mfi);

             // This FAB has cut cells
             amrex::ParallelFor(ubx_grown, [apx_fab,upls_fab,umns_fab,ccvel_fab,xslopes_fab] 
               AMREX_GPU_DEVICE (int i, int j, int k) noexcept
             {
                 // X-faces
                 if (apx_fab(i,j,k) > 0.0)
                 {
                    upls_fab(i,j,k) = ccvel_fab(i  ,j,k,0) - 0.5 * xslopes_fab(i  ,j,k,0);
                    umns_fab(i,j,k) = ccvel_fab(i-1,j,k,0) + 0.5 * xslopes_fab(i-1,j,k,0);
                 }
             });


             amrex::ParallelFor(vbx_grown, [apy_fab,vpls_fab,vmns_fab,ccvel_fab,yslopes_fab] 
               AMREX_GPU_DEVICE (int i, int j, int k) noexcept
             {
                 // Y-faces
                 if (apy_fab(i,j,k) > 0.0)
                 {
                    vpls_fab(i,j,k) = ccvel_fab(i,j  ,k,1) - 0.5 * yslopes_fab(i,j  ,k,1);
                    vmns_fab(i,j,k) = ccvel_fab(i,j-1,k,1) + 0.5 * yslopes_fab(i,j-1,k,1);
                 }
             });

             amrex::ParallelFor(wbx_grown, [apz_fab,wpls_fab,wmns_fab,ccvel_fab,zslopes_fab] 
               AMREX_GPU_DEVICE (int i, int j, int k) noexcept
             {
                 // Z-faces
                 if (apz_fab(i,j,k) > 0.0) {
                    wpls_fab(i,j,k) = ccvel_fab(i,j,k  ,2) - 0.5 * zslopes_fab(i,j,k  ,2);
                    wmns_fab(i,j,k) = ccvel_fab(i,j,k-1,2) + 0.5 * zslopes_fab(i,j,k-1,2);
                 }
             });

          } // Cut cells
#endif
       } // MFIter

       // ****************************************************************************
       // Make sure to fill face-centered values outside our grid before interpolating
       // ****************************************************************************
       upls.FillBoundary(geom[lev].periodicity());
       umns.FillBoundary(geom[lev].periodicity());
       vpls.FillBoundary(geom[lev].periodicity());
       vmns.FillBoundary(geom[lev].periodicity());
       wpls.FillBoundary(geom[lev].periodicity());
       wmns.FillBoundary(geom[lev].periodicity());

#ifdef AMREX_USE_EB
       // ****************************************************************************
       // Do interpolation to centroids -- only for cut cells
       // ****************************************************************************
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
       for (MFIter mfi(*vel_in[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
       {
           // Tilebox
          Box  bx = mfi.tilebox();
          Box ubx = mfi.nodaltilebox(0);
          Box vbx = mfi.nodaltilebox(1);
          Box wbx = mfi.nodaltilebox(2);
      
          // Check efficiently if this tile contains any eb stuff

          const EBFArrayBox&  vel_fab = static_cast<EBFArrayBox const&>((*vel_in[lev])[mfi]);
          const EBCellFlagFab&  flags = vel_fab.getEBCellFlagFab();

          if ( !(flags.getType(amrex::grow(bx,0)) == FabType::covered ||
                 flags.getType(amrex::grow(bx,1)) == FabType::regular ) )
          {
             // Face-centered velocity components
             const auto& umac_fab = (m_u_mac[lev])->array(mfi);
             const auto& vmac_fab = (m_v_mac[lev])->array(mfi);
             const auto& wmac_fab = (m_w_mac[lev])->array(mfi);

             // Face-centered left and right states
             const auto& upls_fab = upls.array(mfi);
             const auto& vpls_fab = vpls.array(mfi);
             const auto& wpls_fab = wpls.array(mfi);
             const auto& umns_fab = umns.array(mfi);
             const auto& vmns_fab = vmns.array(mfi);
             const auto& wmns_fab = wmns.array(mfi);

             // Face-centered areas
             const auto& apx_fab = areafrac[0]->array(mfi);
             const auto& apy_fab = areafrac[1]->array(mfi);
             const auto& apz_fab = areafrac[2]->array(mfi);

             // Face centroids
             const auto& fcx_fab = facecent[0]->array(mfi);
             const auto& fcy_fab = facecent[1]->array(mfi);
             const auto& fcz_fab = facecent[2]->array(mfi);

             const auto& ccm_fab = cc_mask.const_array(mfi);

             amrex::ParallelFor(ubx,[small_vel,huge_vel,apx_fab,ccm_fab,fcx_fab,upls_fab,umns_fab,umac_fab] 
               AMREX_GPU_DEVICE (int i, int j, int k) noexcept
             {
                if (apx_fab(i,j,k) == 0.0)

                    umac_fab(i,j,k) = huge_vel;

                else if (apx_fab(i,j,k) < 1.0) {

                    int jj = j + static_cast<int>(std::copysign(1.0, fcx_fab(i,j,k,0)));
                    int kk = k + static_cast<int>(std::copysign(1.0, fcx_fab(i,j,k,1)));

                    Real fracy = (ccm_fab(i-1,jj,k) || ccm_fab(i,jj,k)) ? std::abs(fcx_fab(i,j,k,0)) : 0.0;
                    Real fracz = (ccm_fab(i-1,j,kk) || ccm_fab(i,j,kk)) ? std::abs(fcx_fab(i,j,k,1)) : 0.0;

                    Real upls_on_centroid = (1.0-fracy)*(1.0-fracz)*upls_fab(i, j,k )+
                                                 fracy *(1.0-fracz)*upls_fab(i,jj,k )+
                                                 fracz *(1.0-fracy)*upls_fab(i, j,kk)+
                                                 fracy *     fracz *upls_fab(i,jj,kk);
                    Real umns_on_centroid = (1.0-fracy)*(1.0-fracz)*umns_fab(i, j,k )+
                                                 fracy *(1.0-fracz)*umns_fab(i,jj,k )+
                                                 fracz *(1.0-fracy)*umns_fab(i, j,kk)+
                                                 fracy *     fracz *umns_fab(i,jj,kk);

                    if ( umns_on_centroid < 0.0 && upls_on_centroid > 0.0 ) {
                       umac_fab(i,j,k) = 0.0;
                    } else {
                       Real avg = 0.5 * ( upls_on_centroid + umns_on_centroid );
                       if ( std::abs(avg) <  small_vel) { umac_fab(i,j,k) = 0.0;
                       } else if (avg >= 0)             { umac_fab(i,j,k) = umns_on_centroid;
                       } else                           { umac_fab(i,j,k) = upls_on_centroid;
                       }
                    }

                 } else {

                    if ( umns_fab(i,j,k) < 0.0 && upls_fab(i,j,k) > 0.0 ) {
                       umac_fab(i,j,k) = 0.0;
                    } else {
                       Real avg = 0.5 * ( upls_fab(i,j,k) + umns_fab(i,j,k) );
                       if ( std::abs(avg) <  small_vel) { umac_fab(i,j,k) = 0.0;
                       } else if (avg >= 0)             { umac_fab(i,j,k) = umns_fab(i,j,k);
                       } else                           { umac_fab(i,j,k) = upls_fab(i,j,k);
                       }
                    }
                 }

             });

             amrex::ParallelFor(vbx,[small_vel,huge_vel,apy_fab,ccm_fab,fcy_fab,vpls_fab,vmns_fab,vmac_fab] 
               AMREX_GPU_DEVICE (int i, int j, int k) noexcept
             {
                if (apy_fab(i,j,k) == 0.0) {

                    vmac_fab(i,j,k) = huge_vel;

                } else if (apy_fab(i,j,k) < 1.0) {

                    int ii = i + static_cast<int>(std::copysign(1.0,fcy_fab(i,j,k,0)));
                    int kk = k + static_cast<int>(std::copysign(1.0,fcy_fab(i,j,k,1)));

                    Real fracx = (ccm_fab(ii,j-1,k ) || ccm_fab(ii,j,k)) ? std::abs(fcy_fab(i,j,k,0)) : 0.0;
                    Real fracz = (ccm_fab( i,j-1,kk) || ccm_fab(i,j,kk)) ? std::abs(fcy_fab(i,j,k,1)) : 0.0;

                    Real vpls_on_centroid = (1.0-fracx)*(1.0-fracz)*vpls_fab(i ,j,k )+
                                                 fracx *(1.0-fracz)*vpls_fab(ii,j,k )+
                                                 fracz *(1.0-fracx)*vpls_fab(i ,j,kk)+
                                                 fracx *     fracz *vpls_fab(ii,j,kk);
                    Real vmns_on_centroid = (1.0-fracx)*(1.0-fracz)*vmns_fab(i ,j,k )+
                                                 fracx *(1.0-fracz)*vmns_fab(ii,j,k )+
                                                 fracz *(1.0-fracx)*vmns_fab(i ,j,kk)+
                                                 fracx *     fracz *vmns_fab(ii,j,kk);

                    if ( vmns_on_centroid < 0.0 && vpls_on_centroid > 0.0 ) {
                       vmac_fab(i,j,k) = 0.0;
                    } else {
                       Real avg = 0.5 * ( vpls_on_centroid + vmns_on_centroid );
                       if ( std::abs(avg) <  small_vel) { vmac_fab(i,j,k) = 0.0;
                       } else if (avg >= 0)             { vmac_fab(i,j,k) = vmns_on_centroid;
                       } else                           { vmac_fab(i,j,k) = vpls_on_centroid;
                       }
                    }

                 } else {

                    if ( vmns_fab(i,j,k) < 0.0 && vpls_fab(i,j,k) > 0.0 ) {
                       vmac_fab(i,j,k) = 0.0;
                    } else {
                       Real avg = 0.5 * ( vpls_fab(i,j,k) + vmns_fab(i,j,k) );
                       if ( std::abs(avg) <  small_vel) { vmac_fab(i,j,k) = 0.0;
                       } else if (avg >= 0)             { vmac_fab(i,j,k) = vmns_fab(i,j,k);
                       } else                           { vmac_fab(i,j,k) = vpls_fab(i,j,k);
                       }
                    }
                 }
             });

             amrex::ParallelFor(wbx,[small_vel,huge_vel,apz_fab,ccm_fab,fcz_fab,wpls_fab,wmns_fab,wmac_fab] 
               AMREX_GPU_DEVICE (int i, int j, int k) noexcept
             {
                if (apz_fab(i,j,k) == 0.0)

                    wmac_fab(i,j,k) = huge_vel;

                else if (apz_fab(i,j,k) < 1.0) {

                    int ii = i + static_cast<int>(std::copysign(1.0,fcz_fab(i,j,k,0)));
                    int jj = j + static_cast<int>(std::copysign(1.0,fcz_fab(i,j,k,1)));

                    Real fracx = (ccm_fab(ii,j,k-1) || ccm_fab(ii,j,k)) ? std::abs(fcz_fab(i,j,k,0)) : 0.0;
                    Real fracy = (ccm_fab(i,jj,k-1) || ccm_fab(i,jj,k)) ? std::abs(fcz_fab(i,j,k,1)) : 0.0;

                    Real wpls_on_centroid = (1.0-fracx)*(1.0-fracy)*wpls_fab(i ,j ,k)+
                                                 fracx *(1.0-fracy)*wpls_fab(ii,j ,k)+
                                                 fracy *(1.0-fracx)*wpls_fab(i ,jj,k)+
                                                 fracx *     fracy *wpls_fab(ii,jj,k);
                    Real wmns_on_centroid = (1.0-fracx)*(1.0-fracy)*wmns_fab(i ,j ,k)+
                                                 fracx *(1.0-fracy)*wmns_fab(ii,j ,k)+
                                                 fracy *(1.0-fracx)*wmns_fab(i ,jj,k)+
                                                 fracx *     fracy *wmns_fab(ii,jj,k);

                    if ( wmns_on_centroid < 0.0 && wpls_on_centroid > 0.0 ) {
                       wmac_fab(i,j,k) = 0.0;
                    } else {
                       Real avg = 0.5 * ( wpls_on_centroid + wmns_on_centroid );
                       if ( std::abs(avg) <  small_vel) { wmac_fab(i,j,k) = 0.0;
                       } else if (avg >= 0)             { wmac_fab(i,j,k) = wmns_on_centroid;
                       } else                           { wmac_fab(i,j,k) = wpls_on_centroid;
                       }
                    }

                 } else {

                    if ( wmns_fab(i,j,k) < 0.0 && wpls_fab(i,j,k) > 0.0 ) {
                       wmac_fab(i,j,k) = 0.0;
                    } else {
                       Real avg = 0.5 * ( wpls_fab(i,j,k) + wmns_fab(i,j,k) );
                       if ( std::abs(avg) <  small_vel) { wmac_fab(i,j,k) = 0.0;
                       } else if (avg >= 0)             { wmac_fab(i,j,k) = wmns_fab(i,j,k);
                       } else                           { wmac_fab(i,j,k) = wpls_fab(i,j,k);
                       }
                    }
                 }
             });

          } // Cut cells
       } // MFIter
#endif

       // Set bcs on u_mac
       set_MAC_velocity_bcs( lev, m_u_mac, m_v_mac, m_w_mac, time );
}
