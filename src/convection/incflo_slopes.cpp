#include <incflo.H>

using namespace amrex;

//
// Compute the slopes of Sborder (velocity, density or tracer)
//
void
incflo::incflo_compute_slopes (int lev, Real time, MultiFab& Sborder,
                               Vector<std::unique_ptr<MultiFab>>& xslopes_in,
                               Vector<std::unique_ptr<MultiFab>>& yslopes_in,
                               Vector<std::unique_ptr<MultiFab>>& zslopes_in,
                               int slopes_comp, int ncomp)
{
    BL_PROFILE("incflo::incflo_compute_slopes");

#ifdef AMREX_USE_EB
    EB_set_covered(Sborder, 0, Sborder.nComp(), 1, covered_val);
#endif

    Box domain(geom[lev].Domain());

    // We initialize slopes to zero in the grown domain ... this is essential
    //    to handle the up-winding at outflow faces
    xslopes_in[lev]->setVal(0.0, slopes_comp, ncomp, xslopes_in[lev]->nGrow());
    yslopes_in[lev]->setVal(0.0, slopes_comp, ncomp, yslopes_in[lev]->nGrow());
    zslopes_in[lev]->setVal(0.0, slopes_comp, ncomp, zslopes_in[lev]->nGrow());

    // ... then set them to this large number in the interior in order to be sure
    //     that no "bad" values go unnoticed
    xslopes_in[lev]->setVal(1.2345e300, slopes_comp, ncomp, 0);
    yslopes_in[lev]->setVal(1.2345e300, slopes_comp, ncomp, 0);
    zslopes_in[lev]->setVal(1.2345e300, slopes_comp, ncomp, 0);

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(Sborder,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
       // Tilebox
       Box bx = mfi.tilebox ();

#ifdef AMREX_USE_EB
       // This is to check efficiently if this tile contains any eb stuff
       const EBFArrayBox&  Sborder_fab = static_cast<EBFArrayBox const&>(Sborder[mfi]);
       const EBCellFlagFab&  flags = Sborder_fab.getEBCellFlagFab();

       if (flags.getType(amrex::grow(bx,0)) != FabType::covered )
#endif
       {
           const auto& state_fab =      Sborder.array(mfi);
           const auto&  xs_fab = xslopes_in[lev]->array(mfi);
           const auto&  ys_fab = yslopes_in[lev]->array(mfi);
           const auto&  zs_fab = zslopes_in[lev]->array(mfi);

#ifdef AMREX_USE_EB
           // No cut cells in tile + 1-cell witdh halo -> use non-eb routine
           if (flags.getType(amrex::grow(bx,1)) == FabType::regular )
#endif
           {
               amrex::ParallelFor(bx,ncomp,
                 [slopes_comp,xs_fab,ys_fab,zs_fab,state_fab] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
               {
                   // X direction
                   Real du_xl = 2.0*(state_fab(i  ,j,k,n) - state_fab(i-1,j,k,n));
                   Real du_xr = 2.0*(state_fab(i+1,j,k,n) - state_fab(i  ,j,k,n));
                   Real du_xc = 0.5*(state_fab(i+1,j,k,n) - state_fab(i-1,j,k,n));

                   Real xslope = amrex::min(std::abs(du_xl),std::abs(du_xc),std::abs(du_xr));
                   xslope          = (du_xr*du_xl > 0.0) ? xslope : 0.0;
                   xs_fab(i,j,k,slopes_comp+n) = (du_xc       > 0.0) ? xslope : -xslope;

                   // Y direction
                   Real du_yl = 2.0*(state_fab(i,j  ,k,n) - state_fab(i,j-1,k,n));
                   Real du_yr = 2.0*(state_fab(i,j+1,k,n) - state_fab(i,j  ,k,n));
                   Real du_yc = 0.5*(state_fab(i,j+1,k,n) - state_fab(i,j-1,k,n));

                   Real yslope = amrex::min(std::abs(du_yl),std::abs(du_yc),std::abs(du_yr));
                   yslope          = (du_yr*du_yl > 0.0) ? yslope : 0.0;
                   ys_fab(i,j,k,slopes_comp+n) = (du_yc       > 0.0) ? yslope : -yslope;

                   // Z direction
                   Real du_zl = 2.0*(state_fab(i,j,k  ,n) - state_fab(i,j,k-1,n));
                   Real du_zr = 2.0*(state_fab(i,j,k+1,n) - state_fab(i,j,k  ,n));
                   Real du_zc = 0.5*(state_fab(i,j,k+1,n) - state_fab(i,j,k-1,n));

                   Real zslope = amrex::min(std::abs(du_zl),std::abs(du_zc),std::abs(du_zr));
                   zslope          = (du_zr*du_zl > 0.0) ? zslope : 0.0;
                   zs_fab(i,j,k,slopes_comp+n) = (du_zc       > 0.0) ? zslope : -zslope;
               });
           }
#ifdef AMREX_USE_EB
           else
           {
               const auto& flag_fab =         flags.array();

               amrex::ParallelFor(bx,ncomp,
                 [slopes_comp,xs_fab,ys_fab,zs_fab,flag_fab,state_fab] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
               {
                   if (flag_fab(i,j,k).isCovered())
                   {
                       xs_fab(i,j,k,slopes_comp+n) = 0.0;
                       ys_fab(i,j,k,slopes_comp+n) = 0.0;
                       zs_fab(i,j,k,slopes_comp+n) = 0.0;
                   }
                   else
                   {
                       // X direction
                       Real du_xl = (flag_fab(i-1,j,k).isCovered()) ? 0.0 :
                           2.0*(state_fab(i  ,j,k,n) - state_fab(i-1,j,k,n));
                       Real du_xr = (flag_fab(i+1,j,k).isCovered()) ? 0.0 :
                           2.0*(state_fab(i+1,j,k,n) - state_fab(i  ,j,k,n));
                       Real du_xc = 0.5*(state_fab(i+1,j,k,n) - state_fab(i-1,j,k,n));

                       Real xslope = amrex::min(std::abs(du_xl),std::abs(du_xc),std::abs(du_xr));
                       xslope          = (du_xr*du_xl > 0.0) ? xslope : 0.0;
                       xs_fab(i,j,k,slopes_comp+n) = (du_xc       > 0.0) ? xslope : -xslope;

                       // Y direction
                       Real du_yl = (flag_fab(i,j-1,k).isCovered()) ? 0.0 :
                           2.0*(state_fab(i,j  ,k,n) - state_fab(i,j-1,k,n));
                       Real du_yr = (flag_fab(i,j+1,k).isCovered()) ? 0.0 :
                           2.0*(state_fab(i,j+1,k,n) - state_fab(i,j  ,k,n));
                       Real du_yc = 0.5*(state_fab(i,j+1,k,n) - state_fab(i,j-1,k,n));

                       Real yslope = amrex::min(std::abs(du_yl),std::abs(du_yc),std::abs(du_yr));
                       yslope          = (du_yr*du_yl > 0.0) ? yslope : 0.0;
                       ys_fab(i,j,k,slopes_comp+n) = (du_yc       > 0.0) ? yslope : -yslope;

                       // Z direction
                       Real du_zl = (flag_fab(i,j,k-1).isCovered()) ? 0.0 :
                           2.0*(state_fab(i,j,k  ,n) - state_fab(i,j,k-1,n));
                       Real du_zr = (flag_fab(i,j,k+1).isCovered()) ? 0.0 :
                           2.0*(state_fab(i,j,k+1,n) - state_fab(i,j,k  ,n));
                       Real du_zc = 0.5*(state_fab(i,j,k+1,n) - state_fab(i,j,k-1,n));

                       Real zslope = amrex::min(std::abs(du_zl),std::abs(du_zc),std::abs(du_zr));
                       zslope          = (du_zr*du_zl > 0.0) ? zslope : 0.0;
                       zs_fab(i,j,k,slopes_comp+n) = (du_zc       > 0.0) ? zslope : -zslope;
                   }
               });
           } // end of cut cell region

           const auto& flag_fab =         flags.array();
#endif
           const int minf = bc_list.get_minf();

           const auto&  ilo_ifab  = bc_ilo[lev]->array();
           const auto&  ihi_ifab  = bc_ihi[lev]->array();
           const auto&  jlo_ifab  = bc_jlo[lev]->array();
           const auto&  jhi_ifab  = bc_jhi[lev]->array();
           const auto&  klo_ifab  = bc_klo[lev]->array();
           const auto&  khi_ifab  = bc_khi[lev]->array();

           amrex::ParallelFor(bx,ncomp,
             [minf,domain,slopes_comp,ilo_ifab,ihi_ifab,jlo_ifab,jhi_ifab,klo_ifab,khi_ifab,
#ifdef AMREX_USE_EB
              xs_fab,ys_fab,zs_fab,flag_fab,state_fab] 
#else
              xs_fab,ys_fab,zs_fab,state_fab] 
#endif
              AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
           {
#ifdef AMREX_USE_EB
               if ( (i == domain.smallEnd(0)) && !flag_fab(i,j,k).isCovered() && ilo_ifab(i-1,j,k,0) == minf)
#else
               if ( (i == domain.smallEnd(0)) &&                                 ilo_ifab(i-1,j,k,0) == minf)
#endif
               {
                   Real du_xl = 2.0*(state_fab(i  ,j,k,n) - state_fab(i-1,j,k,n));
                   Real du_xr = 2.0*(state_fab(i+1,j,k,n) - state_fab(i  ,j,k,n));
                   Real du_xc = (state_fab(i+1,j,k,n)+3.0*state_fab(i,j,k,n)-4.0*state_fab(i-1,j,k,n))/3.0;

                   Real xslope = amrex::min(std::abs(du_xl),std::abs(du_xc),std::abs(du_xr));
                   xslope          = (du_xr*du_xl > 0.0) ? xslope : 0.0;
                   xs_fab(i,j,k,slopes_comp+n) = (du_xc       > 0.0) ? xslope : -xslope;
               }

#ifdef AMREX_USE_EB
               if ( (i == domain.bigEnd(0)) && !flag_fab(i,j,k).isCovered() && ihi_ifab(i+1,j,k,0) == minf)
#else
               if ( (i == domain.bigEnd(0)) &&                                 ihi_ifab(i+1,j,k,0) == minf)
#endif
               {
                   Real du_xl = 2.0*(state_fab(i  ,j,k,n) - state_fab(i-1,j,k,n));
                   Real du_xr = 2.0*(state_fab(i+1,j,k,n) - state_fab(i  ,j,k,n));
                   Real du_xc = -(state_fab(i-1,j,k,n)+3.0*state_fab(i,j,k,n)-4.0*state_fab(i+1,j,k,n))/3.0;

                   Real xslope = amrex::min(std::abs(du_xl),std::abs(du_xc),std::abs(du_xr));
                   xslope          = (du_xr*du_xl > 0.0) ? xslope : 0.0;
                   xs_fab(i,j,k,slopes_comp+n) = (du_xc       > 0.0) ? xslope : -xslope;
               }

#ifdef AMREX_USE_EB
               if ( (j == domain.smallEnd(1)) && !flag_fab(i,j,k).isCovered() && jlo_ifab(i,j-1,k,0) == minf)
#else
               if ( (j == domain.smallEnd(1)) &&                                 jlo_ifab(i,j-1,k,0) == minf)
#endif
               {
                   Real du_yl = 2.0*(state_fab(i,j  ,k,n) - state_fab(i,j-1,k,n));
                   Real du_yr = 2.0*(state_fab(i,j+1,k,n) - state_fab(i,j  ,k,n));
                   Real du_yc = (state_fab(i,j+1,k,n)+3.0*state_fab(i,j,k,n)-4.0*state_fab(i,j-1,k,n))/3.0;

                   Real yslope = amrex::min(std::abs(du_yl),std::abs(du_yc),std::abs(du_yr));
                   yslope          = (du_yr*du_yl > 0.0) ? yslope : 0.0;
                   ys_fab(i,j,k,slopes_comp+n) = (du_yc       > 0.0) ? yslope : -yslope;
               }

#ifdef AMREX_USE_EB
               if ( (j == domain.bigEnd(1)) && !flag_fab(i,j,k).isCovered() && jhi_ifab(i,j+1,k,0) == minf)
#else
               if ( (j == domain.bigEnd(1)) &&                                 jhi_ifab(i,j+1,k,0) == minf)
#endif
               {
                   Real du_yl = 2.0*(state_fab(i,j  ,k,n) - state_fab(i,j-1,k,n));
                   Real du_yr = 2.0*(state_fab(i,j+1,k,n) - state_fab(i,j  ,k,n));
                   Real du_yc = -(state_fab(i,j-1,k,n)+3.0*state_fab(i,j,k,n)-4.0*state_fab(i,j+1,k,n))/3.0;

                   Real yslope = amrex::min(std::abs(du_yl),std::abs(du_yc),std::abs(du_yr));
                   yslope          = (du_yr*du_yl > 0.0) ? yslope : 0.0;
                   ys_fab(i,j,k,slopes_comp+n) = (du_yc       > 0.0) ? yslope : -yslope;
               }

#ifdef AMREX_USE_EB
               if ( (k == domain.smallEnd(2)) && !flag_fab(i,j,k).isCovered() && klo_ifab(i,j,k-1,0) == minf)
#else
               if ( (k == domain.smallEnd(2)) &&                                 klo_ifab(i,j,k-1,0) == minf)
#endif
               {
                   Real du_zl = 2.0*(state_fab(i,j,k  ,n) - state_fab(i,j,k-1,n));
                   Real du_zr = 2.0*(state_fab(i,j,k+1,n) - state_fab(i,j,k  ,n));
                   Real du_zc = (state_fab(i,j,k+1,n)+3.0*state_fab(i,j,k,n)-4.0*state_fab(i,j,k-1,n))/3.0;

                   Real zslope = amrex::min(std::abs(du_zl),std::abs(du_zc),std::abs(du_zr));
                   zslope          = (du_zr*du_zl > 0.0) ? zslope : 0.0;
                   zs_fab(i,j,k,slopes_comp+n) = (du_zc       > 0.0) ? zslope : -zslope;
               }
#ifdef AMREX_USE_EB
               if ( (k == domain.bigEnd(2)) && !flag_fab(i,j,k).isCovered() && khi_ifab(i,j,k+1,0) == minf)
#else
               if ( (k == domain.bigEnd(2)) &&                                 khi_ifab(i,j,k+1,0) == minf)
#endif
               {
                   Real du_zl = 2.0*(state_fab(i,j,k  ,n) - state_fab(i,j,k-1,n));
                   Real du_zr = 2.0*(state_fab(i,j,k+1,n) - state_fab(i,j,k  ,n));
                   Real du_zc = -(state_fab(i,j,k-1,n)+3.0*state_fab(i,j,k,n)-4.0*state_fab(i,j,k+1,n))/3.0;

                   Real zslope = amrex::min(std::abs(du_zl),std::abs(du_zc),std::abs(du_zr));
                   zslope          = (du_zr*du_zl > 0.0) ? zslope : 0.0;
                   zs_fab(i,j,k,slopes_comp+n) = (du_zc       > 0.0) ? zslope : -zslope;
               }
           });
           
        } // not covered
    } // MFIter

    xslopes_in[lev] -> FillBoundary(geom[lev].periodicity());
    yslopes_in[lev] -> FillBoundary(geom[lev].periodicity());
    zslopes_in[lev] -> FillBoundary(geom[lev].periodicity());
}
