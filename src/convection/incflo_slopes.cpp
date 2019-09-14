#include <AMReX_EBMultiFabUtil.H>

#include <incflo.H>
#include <mac_F.H>
#include <convection_F.H>
#include <incflo_divop_conv.hpp>
#include <param_mod_F.H>

//
// Compute the slopes of each velocity component in all three directions.
//
void incflo::ComputeSlopes(int lev, MultiFab& Sborder, 
                           Vector<std::unique_ptr<MultiFab>>& xslopes_in,
                           Vector<std::unique_ptr<MultiFab>>& yslopes_in,
                           Vector<std::unique_ptr<MultiFab>>& zslopes_in,
                           int slopes_comp)
{
    BL_PROFILE("incflo::ComputeVelocitySlopes");

    EB_set_covered(Sborder, covered_val);

    Box domain(geom[lev].Domain());

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for(MFIter mfi(Sborder, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        // Tilebox
        Box bx = mfi.tilebox();

        // this is to check efficiently if this tile contains any eb stuff
        const EBFArrayBox& vel_in_fab = static_cast<EBFArrayBox const&>(Sborder[mfi]);
        const EBCellFlagFab& flags = vel_in_fab.getEBCellFlagFab();

        if(flags.getType(amrex::grow(bx, 0)) == FabType::covered)
        {
		// If tile is completely covered by EB geometry, set slopes
		// value to some very large number so we know if
		// we accidentaly use these covered slopes later in calculations
		xslopes_in[lev]->setVal(1.2345e300, bx, 0, AMREX_SPACEDIM);
		yslopes_in[lev]->setVal(1.2345e300, bx, 0, AMREX_SPACEDIM);
		zslopes_in[lev]->setVal(1.2345e300, bx, 0, AMREX_SPACEDIM);
	}
	else
	{
                const auto&  vel_fab =      Sborder.array(mfi);
                const auto&   xs_fab = xslopes_in[lev]->array(mfi);
                const auto&   ys_fab = yslopes_in[lev]->array(mfi);
                const auto&   zs_fab = zslopes_in[lev]->array(mfi);

                int ncomp = Sborder.nComp();

        	// No cut cells in tile + 1-cell witdh halo -> use non-eb routine
		if(flags.getType(amrex::grow(bx, 1)) == FabType::regular)
		{
                AMREX_HOST_DEVICE_FOR_4D(bx, ncomp, i, j, k, dir,
                {
                   // X direction
                   Real du_xl = 2.0*(vel_fab(i  ,j,k,dir) - vel_fab(i-1,j,k,dir));
                   Real du_xr = 2.0*(vel_fab(i+1,j,k,dir) - vel_fab(i  ,j,k,dir));
                   Real du_xc = 0.5*(vel_fab(i+1,j,k,dir) - vel_fab(i-1,j,k,dir));

                   Real xslope = amrex::min(std::abs(du_xl),std::abs(du_xc),std::abs(du_xr));
                   xslope            = (du_xr*du_xl > 0.0) ? xslope : 0.0;
                   xs_fab(i,j,k,slopes_comp+dir) = (du_xc       > 0.0) ? xslope : -xslope;

                   // Y direction
                   Real du_yl = 2.0*(vel_fab(i,j  ,k,dir) - vel_fab(i,j-1,k,dir));
                   Real du_yr = 2.0*(vel_fab(i,j+1,k,dir) - vel_fab(i,j  ,k,dir));
                   Real du_yc = 0.5*(vel_fab(i,j+1,k,dir) - vel_fab(i,j-1,k,dir));

                   Real yslope = amrex::min(std::abs(du_yl),std::abs(du_yc),std::abs(du_yr));
                   yslope            = (du_yr*du_yl > 0.0) ? yslope : 0.0;
                   ys_fab(i,j,k,slopes_comp+dir) = (du_yc       > 0.0) ? yslope : -yslope;

                   // Z direction
                   Real du_zl = 2.0*(vel_fab(i,j,k  ,dir) - vel_fab(i,j,k-1,dir));
                   Real du_zr = 2.0*(vel_fab(i,j,k+1,dir) - vel_fab(i,j,k  ,dir));
                   Real du_zc = 0.5*(vel_fab(i,j,k+1,dir) - vel_fab(i,j,k-1,dir));

                   Real zslope = amrex::min(std::abs(du_zl),std::abs(du_zc),std::abs(du_zr));
                   zslope            = (du_zr*du_zl > 0.0) ? zslope : 0.0;
                   zs_fab(i,j,k,slopes_comp+dir) = (du_zc       > 0.0) ? zslope : -zslope;
                });

		} else {

                const auto& flag_fab = flags.array();

                AMREX_HOST_DEVICE_FOR_4D(bx, ncomp, i, j, k, dir,
                {
                    if (flag_fab(i,j,k).isCovered())
                    {
                        xs_fab(i,j,k,slopes_comp+dir) = 0.0;
                        ys_fab(i,j,k,slopes_comp+dir) = 0.0;
                        zs_fab(i,j,k,slopes_comp+dir) = 0.0;

                    } else {

                        // X direction
                        Real du_xl = (flag_fab(i-1,j,k).isCovered()) ? 0.0 :
                                     2.0*(vel_fab(i  ,j,k,dir) - vel_fab(i-1,j,k,dir));
                        Real du_xr = (flag_fab(i+1,j,k).isCovered()) ? 0.0 :
                                     2.0*(vel_fab(i+1,j,k,dir) - vel_fab(i  ,j,k,dir));
                        Real du_xc = 0.5*(vel_fab(i+1,j,k,dir) - vel_fab(i-1,j,k,dir));

                        Real xslope = amrex::min(std::abs(du_xl),std::abs(du_xc),std::abs(du_xr));
                        xslope            = (du_xr*du_xl > 0.0) ? xslope : 0.0;
                        xs_fab(i,j,k,slopes_comp+dir) = (du_xc       > 0.0) ? xslope : -xslope;

                        // Y direction
                        Real du_yl = (flag_fab(i,j-1,k).isCovered()) ? 0.0 :
                                     2.0*(vel_fab(i,j  ,k,dir) - vel_fab(i,j-1,k,dir));
                        Real du_yr = (flag_fab(i,j+1,k).isCovered()) ? 0.0 :
                                     2.0*(vel_fab(i,j+1,k,dir) - vel_fab(i,j  ,k,dir));
                        Real du_yc = 0.5*(vel_fab(i,j+1,k,dir) - vel_fab(i,j-1,k,dir));

                        Real yslope = amrex::min(std::abs(du_yl),std::abs(du_yc),std::abs(du_yr));
                        yslope            = (du_yr*du_yl > 0.0) ? yslope : 0.0;
                        ys_fab(i,j,k,slopes_comp+dir) = (du_yc       > 0.0) ? yslope : -yslope;

                        // Z direction
                        Real du_zl = (flag_fab(i,j,k-1).isCovered()) ? 0.0 :
                                     2.0*(vel_fab(i,j,k  ,dir) - vel_fab(i,j,k-1,dir));
                        Real du_zr = (flag_fab(i,j,k+1).isCovered()) ? 0.0 :
                                     2.0*(vel_fab(i,j,k+1,dir) - vel_fab(i,j,k  ,dir));
                        Real du_zc = 0.5*(vel_fab(i,j,k+1,dir) - vel_fab(i,j,k-1,dir));

                        Real zslope = amrex::min(std::abs(du_zl),std::abs(du_zc),std::abs(du_zr));
                        zslope          = (du_zr*du_zl > 0.0) ? zslope : 0.0;
                        zs_fab(i,j,k,slopes_comp+dir) = (du_zc       > 0.0) ? zslope : -zslope;
                    }
                });
	        }

            // TODO -- do we have domain and ilo_fab, etc on GPU???
            // TODO -- we need to use "MINF" from the Fortran, not hard-wire this to 20

            const auto& flag_fab =         flags.array();

            const auto&  ilo_ifab  = bc_ilo[lev]->array();
            const auto&  ihi_ifab  = bc_ihi[lev]->array();
            const auto&  jlo_ifab  = bc_jlo[lev]->array();
            const auto&  jhi_ifab  = bc_jhi[lev]->array();
            const auto&  klo_ifab  = bc_klo[lev]->array();
            const auto&  khi_ifab  = bc_khi[lev]->array();

            AMREX_HOST_DEVICE_FOR_4D(bx, ncomp, i, j, k, dir,
            {
                if ( (i == domain.smallEnd(0)) && !flag_fab(i,j,k).isCovered() && ilo_ifab(i-1,j,k,0) == 20)
                {
                    Real du_xl = 2.0*(vel_fab(i  ,j,k,dir) - vel_fab(i-1,j,k,dir));
                    Real du_xr = 2.0*(vel_fab(i+1,j,k,dir) - vel_fab(i  ,j,k,dir));
                    Real du_xc = (vel_fab(i+1,j,k,dir)+3.0*vel_fab(i,j,k,dir)-4.0*vel_fab(i-1,j,k,dir))/3.0;

                    Real xslope = amrex::min(std::abs(du_xl),std::abs(du_xc),std::abs(du_xr));
                    xslope            = (du_xr*du_xl > 0.0) ? xslope : 0.0;
                    xs_fab(i,j,k,slopes_comp+dir) = (du_xc       > 0.0) ? xslope : -xslope;
                }
                if ( (i == domain.bigEnd(0)) && !flag_fab(i,j,k).isCovered() && ihi_ifab(i+1,j,k,0) == 20)
                {
                    Real du_xl = 2.0*(vel_fab(i  ,j,k,dir) - vel_fab(i-1,j,k,dir));
                    Real du_xr = 2.0*(vel_fab(i+1,j,k,dir) - vel_fab(i  ,j,k,dir));
                    Real du_xc = -(vel_fab(i-1,j,k,dir)+3.0*vel_fab(i,j,k,dir)-4.0*vel_fab(i+1,j,k,dir))/3.0;

                    Real xslope = amrex::min(std::abs(du_xl),std::abs(du_xc),std::abs(du_xr));
                    xslope            = (du_xr*du_xl > 0.0) ? xslope : 0.0;
                    xs_fab(i,j,k,slopes_comp+dir) = (du_xc       > 0.0) ? xslope : -xslope;
                }

                if ( (j == domain.smallEnd(1)) && !flag_fab(i,j,k).isCovered() && jlo_ifab(i,j-1,k,0) == 20)
                {
                    Real du_yl = 2.0*(vel_fab(i,j  ,k,dir) - vel_fab(i,j-1,k,dir));
                    Real du_yr = 2.0*(vel_fab(i,j+1,k,dir) - vel_fab(i,j  ,k,dir));
                    Real du_yc = (vel_fab(i,j+1,k,dir)+3.0*vel_fab(i,j,k,dir)-4.0*vel_fab(i,j-1,k,dir))/3.0;

                    Real yslope = amrex::min(std::abs(du_yl),std::abs(du_yc),std::abs(du_yr));
                    yslope            = (du_yr*du_yl > 0.0) ? yslope : 0.0;
                    ys_fab(i,j,k,slopes_comp+dir) = (du_yc       > 0.0) ? yslope : -yslope;
                }
                if ( (j == domain.bigEnd(1)) && !flag_fab(i,j,k).isCovered() && jhi_ifab(i,j+1,k,0) == 20)
                {
                    Real du_yl = 2.0*(vel_fab(i,j  ,k,dir) - vel_fab(i,j-1,k,dir));
                    Real du_yr = 2.0*(vel_fab(i,j+1,k,dir) - vel_fab(i,j  ,k,dir));
                    Real du_yc = -(vel_fab(i,j-1,k,dir)+3.0*vel_fab(i,j,k,dir)-4.0*vel_fab(i,j+1,k,dir))/3.0;

                    Real yslope = amrex::min(std::abs(du_yl),std::abs(du_yc),std::abs(du_yr));
                    yslope            = (du_yr*du_yl > 0.0) ? yslope : 0.0;
                    ys_fab(i,j,k,slopes_comp+dir) = (du_yc       > 0.0) ? yslope : -yslope;
                }

                if ( (k == domain.smallEnd(2)) && !flag_fab(i,j,k).isCovered() && klo_ifab(i,j,k-1,0) == 20)
                {
                    Real du_zl = 2.0*(vel_fab(i,j,k  ,dir) - vel_fab(i,j,k-1,dir));
                    Real du_zr = 2.0*(vel_fab(i,j,k+1,dir) - vel_fab(i,j,k  ,dir));
                    Real du_zc = (vel_fab(i,j,k+1,dir)+3.0*vel_fab(i,j,k,dir)-4.0*vel_fab(i,j,k-1,dir))/3.0;

                    Real zslope = amrex::min(std::abs(du_zl),std::abs(du_zc),std::abs(du_zr));
                    zslope            = (du_zr*du_zl > 0.0) ? zslope : 0.0;
                    zs_fab(i,j,k,slopes_comp+dir) = (du_zc       > 0.0) ? zslope : -zslope;
                }
                if ( (k == domain.bigEnd(2)) && !flag_fab(i,j,k).isCovered() && khi_ifab(i,j,k+1,0) == 20)
                {
                    Real du_zl = 2.0*(vel_fab(i,j,k  ,dir) - vel_fab(i,j,k-1,dir));
                    Real du_zr = 2.0*(vel_fab(i,j,k+1,dir) - vel_fab(i,j,k  ,dir));
                    Real du_zc = -(vel_fab(i,j,k-1,dir)+3.0*vel_fab(i,j,k,dir)-4.0*vel_fab(i,j,k+1,dir))/3.0;

                    Real zslope = amrex::min(std::abs(du_zl),std::abs(du_zc),std::abs(du_zr));
                    zslope            = (du_zr*du_zl > 0.0) ? zslope : 0.0;
                    zs_fab(i,j,k,slopes_comp+dir) = (du_zc       > 0.0) ? zslope : -zslope;
                }
            });
	}
    }

    xslopes_in[lev]->FillBoundary(geom[lev].periodicity());
    yslopes_in[lev]->FillBoundary(geom[lev].periodicity());
    zslopes_in[lev]->FillBoundary(geom[lev].periodicity());
}
