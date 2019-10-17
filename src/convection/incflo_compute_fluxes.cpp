#include <incflo.H>
#include <param_mod_F.H>

#include <AMReX_REAL.H>
#include <AMReX_BLFort.H>
#include <AMReX_SPACE.H>
#include <AMReX_Array.H>

namespace ugradu_aux {

//
// Compute upwind non-normal velocity
//
AMREX_GPU_HOST_DEVICE
Real
upwind(const Real velocity_minus,
       const Real velocity_plus,
       const Real u_edge)
{
  // Small value to protect against tiny velocities used in upwinding
  const Real small_velocity(1.e-10);

  if(std::abs(u_edge) < small_velocity)
    return .5*(velocity_minus+velocity_plus);

  return u_edge > 0 ? velocity_minus : velocity_plus;
}

AMREX_GPU_HOST_DEVICE
bool
is_equal_to_any(const int bc,
                const int* bc_types,
                const int size)
{
  for(int i(0); i < size; ++i)
  {
    if(bc == bc_types[i])
      return true;
  }
  return false;
}

} // end namespace ugradu_aux

using namespace ugradu_aux;

//
// Compute the three components of the convection term
//
void
incflo::incflo_compute_fluxes(int lev,
                          Vector< std::unique_ptr<MultiFab> >& a_fx,
                          Vector< std::unique_ptr<MultiFab> >& a_fy,
                          Vector< std::unique_ptr<MultiFab> >& a_fz,
                          Vector< std::unique_ptr<MultiFab> >& state_in,
                          const int state_comp, const int ncomp,
                          Vector< std::unique_ptr<MultiFab> >& xslopes_in,
                          Vector< std::unique_ptr<MultiFab> >& yslopes_in,
                          Vector< std::unique_ptr<MultiFab> >& zslopes_in,
                          const int slopes_comp,
                          Vector< std::unique_ptr<MultiFab> >& u_mac,
                          Vector< std::unique_ptr<MultiFab> >& v_mac,
                          Vector< std::unique_ptr<MultiFab> >& w_mac)
{
        Box domain(geom[lev].Domain());

        // Get EB geometric info
        Array< const MultiCutFab*,AMREX_SPACEDIM> areafrac;
        Array< const MultiCutFab*,AMREX_SPACEDIM> facecent;
        const amrex::MultiFab*                    volfrac;
        const amrex::MultiCutFab*                 bndrycent;

        areafrac  =   ebfactory[lev] -> getAreaFrac();
        facecent  =   ebfactory[lev] -> getFaceCent();
        volfrac   = &(ebfactory[lev] -> getVolFrac());
        bndrycent = &(ebfactory[lev] -> getBndryCent());

        // Create cc_mask
        iMultiFab cc_mask(grids[lev], dmap[lev], 1, 1);

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
       {
           std::vector< std::pair<int,Box> > isects;
           const std::vector<IntVect>& pshifts = geom[lev].periodicity().shiftIntVect();
           const BoxArray& ba = cc_mask.boxArray();
           for (MFIter mfi(cc_mask); mfi.isValid(); ++mfi)
           {
               Array4<int> const& fab = cc_mask.array(mfi);

               const Box& bx = mfi.fabbox();
               for (const auto& iv : pshifts)
               {
                   ba.intersections(bx+iv, isects);
                   for (const auto& is : isects)
                   {
                       const Box& b = is.second-iv;
                       AMREX_FOR_3D ( b, i, j, k,
                       {
                           fab(i,j,k) = 1;
                       });
                   }
               }
               // NOTE: here we do not need host-device synchronization since it
               // is already included in the MFIter destructor
           }
        }

        for (MFIter mfi(*state_in[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            // Tilebox
            Box bx = mfi.tilebox ();

            // this is to check efficiently if this tile contains any eb stuff
            const EBFArrayBox&  state_fab = static_cast<EBFArrayBox const&>((*state_in[lev])[mfi]);
            const EBCellFlagFab&  flags = state_fab.getEBCellFlagFab();

            if (flags.getType(amrex::grow(bx,0)) == FabType::covered )
            {
                 const Box ubx = amrex::surroundingNodes(bx,0);
                 const Box vbx = amrex::surroundingNodes(bx,1);
                 const Box wbx = amrex::surroundingNodes(bx,2);
                 a_fx[lev]->setVal(covered_val, ubx, 0, ncomp);
                 a_fy[lev]->setVal(covered_val, vbx, 0, ncomp);
                 a_fz[lev]->setVal(covered_val, wbx, 0, ncomp);
            }
            else
            {
                // No cut cells in tile + nghost-cell witdh halo -> use non-eb routine
                if (flags.getType(amrex::grow(bx,nghost)) == FabType::regular )
                {
                    incflo_compute_ugradu(lev, bx, a_fx, a_fy, a_fz, state_in, state_comp, ncomp,
                                        xslopes_in, yslopes_in, zslopes_in, slopes_comp,
                                        u_mac, v_mac, w_mac, &mfi, domain);
                }
                else
                {
                    incflo_compute_ugradu_eb(lev, bx, a_fx, a_fy, a_fz, state_in, state_comp, ncomp,
                                           xslopes_in, yslopes_in, zslopes_in, slopes_comp,
                                           u_mac, v_mac, w_mac, &mfi, areafrac, facecent,
                                           volfrac, bndrycent, &cc_mask, domain, flags);
                }
            }
        } // MFIter
}

void
incflo::incflo_compute_ugradu( const int lev, Box& bx,
                           Vector< std::unique_ptr<MultiFab> >& a_fx,
                           Vector< std::unique_ptr<MultiFab> >& a_fy,
                           Vector< std::unique_ptr<MultiFab> >& a_fz,
                           Vector< std::unique_ptr<MultiFab> >& state_in,
                           const int state_comp, const int ncomp,
                           Vector< std::unique_ptr<MultiFab> >& xslopes_in,
                           Vector< std::unique_ptr<MultiFab> >& yslopes_in,
                           Vector< std::unique_ptr<MultiFab> >& zslopes_in,
                           const int slopes_comp,
                           Vector< std::unique_ptr<MultiFab> >& u_mac,
                           Vector< std::unique_ptr<MultiFab> >& v_mac,
                           Vector< std::unique_ptr<MultiFab> >& w_mac,
                           MFIter* mfi, Box& domain)
{
  const amrex::Dim3 dom_low = amrex::lbound(domain);
  const amrex::Dim3 dom_high = amrex::ubound(domain);

  Array4<Real> const& state     = state_in[lev]->array(*mfi);

  Array4<Real> const& u = u_mac[lev]->array(*mfi);
  Array4<Real> const& v = v_mac[lev]->array(*mfi);
  Array4<Real> const& w = w_mac[lev]->array(*mfi);

  Array4<Real> const& fx = a_fx[lev]->array(*mfi);
  Array4<Real> const& fy = a_fy[lev]->array(*mfi);
  Array4<Real> const& fz = a_fz[lev]->array(*mfi);

  Array4<Real> const& x_slopes = xslopes_in[lev]->array(*mfi);
  Array4<Real> const& y_slopes = yslopes_in[lev]->array(*mfi);
  Array4<Real> const& z_slopes = zslopes_in[lev]->array(*mfi);

  Array4<int> const& bct_ilo = bc_ilo[lev]->array();
  Array4<int> const& bct_ihi = bc_ihi[lev]->array();
  Array4<int> const& bct_jlo = bc_jlo[lev]->array();
  Array4<int> const& bct_jhi = bc_jhi[lev]->array();
  Array4<int> const& bct_klo = bc_klo[lev]->array();
  Array4<int> const& bct_khi = bc_khi[lev]->array();

  const Box ubx       = amrex::surroundingNodes(bx,0);
  const Box vbx       = amrex::surroundingNodes(bx,1);
  const Box wbx       = amrex::surroundingNodes(bx,2);

  // Vectorize the boundary conditions list in order to use it in lambda
  // functions
  const GpuArray<int, 3> bc_types =
    {bc_list.get_minf(), bc_list.get_pinf(), bc_list.get_pout()};

  AMREX_FOR_4D(ubx, ncomp, i, j, k, n,
  {
    Real state_w(0)  ;
    //Real state_e(0); DECLARED_BUT_NEVER_REFERENCED
    Real state_mns(0); Real state_pls(0);

    //
    // West face
    //
    // In the case of MINF       we are using the prescribed Dirichlet value
    // In the case of PINF, POUT we are using the upwind value
    if ((i == dom_low.x) and
     ugradu_aux::is_equal_to_any(bct_ilo(dom_low.x-1,j,k,0),
                                 bc_types.data(), bc_types.size()))
    {
      state_w = state(i-1,j,k,state_comp+n);

    } else if ((i == dom_high.x+1) and
     ugradu_aux::is_equal_to_any(bct_ihi(dom_high.x+1,j,k,0),
                                 bc_types.data(), bc_types.size()))
    {
      state_w = state(i,j,k,state_comp+n);
    } else {
      state_pls = state(i  ,j,k,state_comp+n) - .5*x_slopes(i  ,j,k,slopes_comp+n);
      state_mns = state(i-1,j,k,state_comp+n) + .5*x_slopes(i-1,j,k,slopes_comp+n);
      state_w = upwind( state_mns, state_pls, u(i,j,k) );
    }
    fx(i,j,k,n) = u(i,j,k) * state_w;
  });

  AMREX_FOR_4D(vbx, ncomp, i, j, k, n,
  {
    Real state_s(0)  ;
    //Real state_n(0); DECLARED_BUT_NEVER_REFERENCED
    Real state_mns(0); Real state_pls(0);

    //
    // South face
    //
    // In the case of MINF       we are using the prescribed Dirichlet value
    // In the case of PINF, POUT we are using the upwind value
    if((j == dom_low.y) and
     ugradu_aux::is_equal_to_any(bct_jlo(i,dom_low.y-1,k,0),
                                 bc_types.data(), bc_types.size()))
    {
      state_s = state(i,j-1,k,state_comp+n);
    } else if ((j == dom_high.y+1) and
     ugradu_aux::is_equal_to_any(bct_jhi(i,dom_high.y+1,k,0),
                                 bc_types.data(), bc_types.size()))
    {
      state_s = state(i,j,k,state_comp+n);
    } else {
      state_pls = state(i,j  ,k,state_comp+n) - .5*y_slopes(i,j  ,k,slopes_comp+n);
      state_mns = state(i,j-1,k,state_comp+n) + .5*y_slopes(i,j-1,k,slopes_comp+n);

      state_s = upwind( state_mns, state_pls, v(i,j,k) );
    }
    fy(i,j,k,n) = v(i,j,k) * state_s;
  });

  AMREX_FOR_4D(wbx, ncomp, i, j, k, n,
  {
    Real state_b(0)  ;
    // Real state_t(0); DECLARED_BUT_NEVER_REFERENCED
    Real state_mns(0); Real state_pls(0);
    //
    // Bottom face
    //
    // In the case of MINF       we are using the prescribed Dirichlet value
    // In the case of PINF, POUT we are using the upwind value
    if((k == dom_low.z) and
     ugradu_aux::is_equal_to_any(bct_klo(i,j,dom_low.z-1,0),
                                 bc_types.data(), bc_types.size()))
    {
      state_b = state(i,j,k-1,state_comp+n);
    } else if ((k == dom_high.z+1) and
     ugradu_aux::is_equal_to_any(bct_khi(i,j,dom_high.z+1,0),
                                 bc_types.data(), bc_types.size()))
    {
      state_b = state(i,j,k,state_comp+n);
    } else {
      state_pls = state(i,j,k  ,state_comp+n) - .5*z_slopes(i,j,k  ,slopes_comp+n);
      state_mns = state(i,j,k-1,state_comp+n) + .5*z_slopes(i,j,k-1,slopes_comp+n);

      state_b = upwind( state_mns, state_pls, w(i,j,k) );
    }
    fz(i,j,k,n) = w(i,j,k) * state_b;
  });

  Gpu::synchronize();
}


//
// Compute the three components of the convection term when we have embedded
// boundaries
//
void
incflo::incflo_compute_ugradu_eb(const int lev, Box& bx,
                             Vector< std::unique_ptr<MultiFab> >& a_fx,
                             Vector< std::unique_ptr<MultiFab> >& a_fy,
                             Vector< std::unique_ptr<MultiFab> >& a_fz,
                             Vector< std::unique_ptr<MultiFab> >& state_in,
                             const int state_comp, const int ncomp,
                             Vector< std::unique_ptr<MultiFab> >& xslopes_in,
                             Vector< std::unique_ptr<MultiFab> >& yslopes_in,
                             Vector< std::unique_ptr<MultiFab> >& zslopes_in,
                             const int slopes_comp,
                             Vector< std::unique_ptr<MultiFab> >& u_mac,
                             Vector< std::unique_ptr<MultiFab> >& v_mac,
                             Vector< std::unique_ptr<MultiFab> >& w_mac,
                             MFIter* mfi,
                             Array<const MultiCutFab*,AMREX_SPACEDIM>& areafrac,
                             Array<const MultiCutFab*,AMREX_SPACEDIM>& facecent,
                             const MultiFab* volfrac,
                             const MultiCutFab* bndrycent,
                             const iMultiFab* cc_mask,
                             Box& domain, const EBCellFlagFab& flags)
{
  const amrex::Dim3 dom_low = amrex::lbound(domain);
  const amrex::Dim3 dom_high = amrex::ubound(domain);

  Array4<Real> const& state     = state_in[lev]->array(*mfi);

  Array4<const Real> const& areafrac_x = areafrac[0]->array(*mfi);
  Array4<const Real> const& areafrac_y = areafrac[1]->array(*mfi);
  Array4<const Real> const& areafrac_z = areafrac[2]->array(*mfi);

  Array4<Real> const& u = u_mac[lev]->array(*mfi);
  Array4<Real> const& v = v_mac[lev]->array(*mfi);
  Array4<Real> const& w = w_mac[lev]->array(*mfi);

  Array4<Real> const& fx = a_fx[lev]->array(*mfi);
  Array4<Real> const& fy = a_fy[lev]->array(*mfi);
  Array4<Real> const& fz = a_fz[lev]->array(*mfi);

  Array4<Real> const& x_slopes = xslopes_in[lev]->array(*mfi);
  Array4<Real> const& y_slopes = yslopes_in[lev]->array(*mfi);
  Array4<Real> const& z_slopes = zslopes_in[lev]->array(*mfi);

  Array4<int> const& bc_ilo_type = bc_ilo[lev]->array();
  Array4<int> const& bc_ihi_type = bc_ihi[lev]->array();
  Array4<int> const& bc_jlo_type = bc_jlo[lev]->array();
  Array4<int> const& bc_jhi_type = bc_jhi[lev]->array();
  Array4<int> const& bc_klo_type = bc_klo[lev]->array();
  Array4<int> const& bc_khi_type = bc_khi[lev]->array();

  const Box ubx       = amrex::surroundingNodes(bx,0);
  const Box vbx       = amrex::surroundingNodes(bx,1);
  const Box wbx       = amrex::surroundingNodes(bx,2);

  const Box ubx_grown = amrex::surroundingNodes(amrex::grow(bx,1),0);
  const Box vbx_grown = amrex::surroundingNodes(amrex::grow(bx,1),1);
  const Box wbx_grown = amrex::surroundingNodes(amrex::grow(bx,1),2);

  FArrayBox s_on_x_face(ubx_grown, ncomp);
  FArrayBox s_on_y_face(vbx_grown, ncomp);
  FArrayBox s_on_z_face(wbx_grown, ncomp);

  Array4<Real> const& sx = s_on_x_face.array();
  Array4<Real> const& sy = s_on_y_face.array();
  Array4<Real> const& sz = s_on_z_face.array();

  // Face centroids
  const auto& fcx_fab = facecent[0]->array(*mfi);
  const auto& fcy_fab = facecent[1]->array(*mfi);
  const auto& fcz_fab = facecent[2]->array(*mfi);

  const auto& ccm_fab = cc_mask->const_array(*mfi);

  const GpuArray<int, 3> bc_types =
    {bc_list.get_minf(), bc_list.get_pinf(), bc_list.get_pout()};

  const Real my_huge = get_my_huge();
  //
  // First compute the convective fluxes at the face center
  // Do this on ALL faces on the tile, i.e. INCLUDE as many ghost faces as
  // possible
  //

  //
  // ===================== X =====================
  //
  AMREX_FOR_4D(ubx_grown, ncomp, i, j, k, n,
  {
    Real upls(0); Real umns(0);

    if( areafrac_x(i,j,k) > 0 ) {
      if( i <= dom_low.x and
       ugradu_aux::is_equal_to_any(bc_ilo_type(dom_low.x-1,j,k,0),
                                   bc_types.data(), bc_types.size()))
      {
        sx(i,j,k,n) = state(dom_low.x-1,j,k,state_comp+n);
      }
      else if( i >= dom_high.x+1 and
       ugradu_aux::is_equal_to_any(bc_ihi_type(dom_high.x+1,j,k,0),
                                   bc_types.data(), bc_types.size()))
      {
        sx(i,j,k,n) = state(dom_high.x+1,j,k,state_comp+n);
      }
      else {
        upls = state(i  ,j,k,state_comp+n) - .5*x_slopes(i  ,j,k,slopes_comp+n);
        umns = state(i-1,j,k,state_comp+n) + .5*x_slopes(i-1,j,k,slopes_comp+n);

        sx(i,j,k,n) = upwind( umns, upls, u(i,j,k) );
      }
    } else {
        sx(i,j,k,n) = my_huge;
    }
  });

  AMREX_FOR_4D(ubx, ncomp, i, j, k, n,
  {
    if( areafrac_x(i,j,k) > 0 ) {
       int jj = j + static_cast<int>(std::copysign(1.0, fcx_fab(i,j,k,0)));
       int kk = k + static_cast<int>(std::copysign(1.0, fcx_fab(i,j,k,1)));

       Real fracy = (ccm_fab(i-1,jj,k) || ccm_fab(i,jj,k)) ? std::abs(fcx_fab(i,j,k,0)) : 0.0;
       Real fracz = (ccm_fab(i-1,j,kk) || ccm_fab(i,j,kk)) ? std::abs(fcx_fab(i,j,k,1)) : 0.0;

       Real s_on_x_centroid = (1.0-fracy)*(1.0-fracz)*sx(i, j,k ,n)+
                                   fracy *(1.0-fracz)*sx(i,jj,k ,n)+
                                   fracz *(1.0-fracy)*sx(i, j,kk,n)+
                                   fracy *     fracz *sx(i,jj,kk,n);

       fx(i,j,k,n) = u(i,j,k) * s_on_x_centroid;
    } else
       fx(i,j,k,n) = my_huge;
  });

  //
  // ===================== Y =====================
  //
  AMREX_FOR_4D(vbx_grown, ncomp, i, j, k, n,
  {
    Real vpls(0); Real vmns(0);

    if( areafrac_y(i,j,k) > 0 ) {
      if( j <= dom_low.y and
       ugradu_aux::is_equal_to_any(bc_jlo_type(i,dom_low.y-1,k,0),
                                   bc_types.data(), bc_types.size()))
      {
        sy(i,j,k,n) = state(i,dom_low.y-1,k,state_comp+n);
      }
      else if( j >= dom_high.y+1 and
       ugradu_aux::is_equal_to_any(bc_jhi_type(i,dom_high.y+1,k,0),
                                   bc_types.data(), bc_types.size()))
      {
        sy(i,j,k,n) = state(i,dom_high.y+1,k,state_comp+n);
      }
      else {
        vpls = state(i,j  ,k,state_comp+n) - .5*y_slopes(i,j  ,k,slopes_comp+n);
        vmns = state(i,j-1,k,state_comp+n) + .5*y_slopes(i,j-1,k,slopes_comp+n);

        sy(i,j,k,n) = upwind( vmns, vpls, v(i,j,k) );
      }
    }
    else {
        sy(i,j,k,n) = my_huge;
    }
  });

  AMREX_FOR_4D(vbx, ncomp, i, j, k, n,
  {
    if ( areafrac_y(i,j,k) > 0 ) {
       int ii = i + static_cast<int>(std::copysign(1.0,fcy_fab(i,j,k,0)));
       int kk = k + static_cast<int>(std::copysign(1.0,fcy_fab(i,j,k,1)));

       Real fracx = (ccm_fab(ii,j-1,k) || ccm_fab(ii,j,k)) ? std::abs(fcy_fab(i,j,k,0)) : 0.0;
       Real fracz = (ccm_fab(i,j-1,kk) || ccm_fab(i,j,kk)) ? std::abs(fcy_fab(i,j,k,1)) : 0.0;

       Real s_on_y_centroid = (1.0-fracx)*(1.0-fracz)*sy(i ,j,k ,n)+
                                   fracx *(1.0-fracz)*sy(ii,j,k ,n)+
                                   fracz *(1.0-fracx)*sy(i ,j,kk,n)+
                                   fracx *     fracz *sy(ii,j,kk,n);
       fy(i,j,k,n) = v(i,j,k) * s_on_y_centroid;
    } else
       fy(i,j,k,n) = my_huge;

  });

  //
  // ===================== Z =====================
  //
  AMREX_FOR_4D(wbx_grown, ncomp, i, j, k, n,
  {
    Real wpls(0); Real wmns(0);

    if( areafrac_z(i,j,k) > 0 ) {
      if( k <= dom_low.z and
       ugradu_aux::is_equal_to_any(bc_klo_type(i,j,dom_low.z-1,0),
                                   bc_types.data(), bc_types.size()))
      {
        sz(i,j,k,n) = state(i,j,dom_low.z-1,state_comp+n);
      }
      else if( k >= dom_high.z+1 and
       ugradu_aux::is_equal_to_any(bc_khi_type(i,j,dom_high.z+1,0),
                                   bc_types.data(), bc_types.size()))
      {
        sz(i,j,k,n) = state(i,j,dom_high.z+1,state_comp+n);
      }
      else {
        wpls = state(i,j,k  ,state_comp+n) - .5*z_slopes(i,j,k  ,slopes_comp+n);
        wmns = state(i,j,k-1,state_comp+n) + .5*z_slopes(i,j,k-1,slopes_comp+n);

        sz(i,j,k,n) = upwind( wmns, wpls, w(i,j,k) );
      }
    }
    else {
        sz(i,j,k,n) = my_huge;
    }
  });

  AMREX_FOR_4D(wbx, ncomp, i, j, k, n,
  {
    if( areafrac_z(i,j,k) > 0 ) {
       int ii = i + static_cast<int>(std::copysign(1.0,fcz_fab(i,j,k,0)));
       int jj = j + static_cast<int>(std::copysign(1.0,fcz_fab(i,j,k,1)));

       Real fracx = (ccm_fab(ii,j,k-1) || ccm_fab(ii,j,k)) ? std::abs(fcz_fab(i,j,k,0)) : 0.0;
       Real fracy = (ccm_fab(i,jj,k-1) || ccm_fab(i,jj,k)) ? std::abs(fcz_fab(i,j,k,1)) : 0.0;

       Real s_on_z_centroid = (1.0-fracx)*(1.0-fracy)*sz(i ,j ,k,n)+
                                   fracx *(1.0-fracy)*sz(ii,j ,k,n)+
                                   fracy *(1.0-fracx)*sz(i ,jj,k,n)+
                                   fracx *     fracy *sz(ii,jj,k,n);

       fz(i,j,k,n) = w(i,j,k) * s_on_z_centroid;
    } else
       fz(i,j,k,n) = my_huge;
  });

  Gpu::synchronize();
}
