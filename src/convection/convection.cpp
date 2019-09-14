#include <AMReX_EBMultiFabUtil.H>

#include <incflo.H>
#include <mac_F.H>
#include <convection_F.H>
#include <incflo_divop_conv.hpp>
#include <param_mod_F.H>

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
// Compute acc using the vel passed in
//
void incflo::ComputeUGradU(Vector<std::unique_ptr<MultiFab>>& conv_in,
                           Vector<std::unique_ptr<MultiFab>>& vel_in,
                           Vector<std::unique_ptr<MultiFab>>& tracer_in,
                           Real time)
{
    BL_PROFILE("incflo::ComputeUGradU");

    // Extrapolate velocity field to cell faces
    ComputeVelocityAtFaces(vel_in, time);

    // Do projection on all AMR-level_ins in one shot
    mac_projection->apply_projection(m_u_mac, m_v_mac, m_w_mac, ro, time, steady_state);

    for(int lev = 0; lev <= finest_level; lev++)
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

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for(MFIter mfi(*vel_in[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            // Tilebox
            Box bx = mfi.tilebox();

            // this is to check efficiently if this tile contains any eb stuff
            const EBFArrayBox& vel_in_fab = static_cast<EBFArrayBox const&>((*vel_in[lev])[mfi]);
            const EBCellFlagFab& flags = vel_in_fab.getEBCellFlagFab();

            if(flags.getType(amrex::grow(bx, 0)) == FabType::covered)
            {
                // If tile is completely covered by EB geometry, set slopes
                // value to some very large number so we know if
                // we accidentaly use these covered slopes later in calculations
                (*conv_in[lev])[mfi].setVal(1.2345e300, bx, 0, AMREX_SPACEDIM);
            }
            else
            {
                // No cut cells in tile + nghost-cell witdh halo -> use non-eb routine
                if(flags.getType(amrex::grow(bx, nghost)) == FabType::regular)
                {
                    incflo_compute_ugradu(bx, conv_in, vel_in, m_u_mac, m_v_mac, m_w_mac, &mfi, domain, lev);

//                  if (advance_tracer)
//                     incflo_compute_ugradu(bx, conv_in, tracer_in, m_u_mac, m_v_mac, m_w_mac, &mfi, domain, lev);
                }
                else
                {

                    incflo_compute_ugradu_eb(bx, conv_in, vel_in, m_u_mac, m_v_mac, m_w_mac, &mfi, areafrac, facecent,
                                             volfrac, bndrycent, domain, flags, lev);

//                  if (advance_tracer)
//                     incflo_compute_ugradu_eb(bx, conv_in, tracer_in, m_u_mac, m_v_mac, m_w_mac, &mfi, areafrac, facecent,
//                                              volfrac, bndrycent, domain, flags, lev);

                }
            }
        }
    }
}

//
// TODO: Documentation
//
void incflo::ComputeVelocityAtFaces(Vector<std::unique_ptr<MultiFab>>& vel_in, Real time)
{
    BL_PROFILE("incflo::ComputeVelocityAtFaces");

    for(int lev = 0; lev <= finest_level; lev++)
    {
        Box domain(geom[lev].Domain());

        // State with ghost cells
        MultiFab Sborder(grids[lev], dmap[lev], vel[lev]->nComp(), nghost,
                         MFInfo(), *ebfactory[lev]);
        FillPatchVel(lev, time, Sborder, 0, Sborder.nComp());

        // Compute the slopes
        ComputeVelocitySlopes(lev, Sborder);

        // Copy each FAB back from Sborder into the vel array, complete with filled ghost cells
        MultiFab::Copy (*vel_in[lev], Sborder, 0, 0, vel_in[lev]->nComp(), vel_in[lev]->nGrow());

        // Get EB geometric info
        Array<const MultiCutFab*, AMREX_SPACEDIM> areafrac;
        areafrac = ebfactory[lev]->getAreaFrac();

    // Then compute velocity at faces
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for(MFIter mfi(Sborder, TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            // Tilebox
            Box bx = mfi.tilebox();
            Box ubx = mfi.tilebox(e_x);
            Box vbx = mfi.tilebox(e_y);
            Box wbx = mfi.tilebox(e_z);

            // this is to check efficiently if this tile contains any eb stuff
            const EBFArrayBox& vel_in_fab = static_cast<EBFArrayBox const&>(Sborder[mfi]);
            const EBCellFlagFab& flags = vel_in_fab.getEBCellFlagFab();

            Real small_vel = 1.e-10;
            Real  huge_vel = 1.e100;

            // Cell-centered velocity
            const auto& ccvel_fab = vel[lev]->array(mfi);

            // Cell-centered slopes
            const auto& xslopes_fab = (xslopes[lev])->array(mfi);
            const auto& yslopes_fab = (yslopes[lev])->array(mfi);
            const auto& zslopes_fab = (zslopes[lev])->array(mfi);

            // Face-centered velocity components
            const auto& umac_fab = (m_u_mac[lev])->array(mfi);
            const auto& vmac_fab = (m_v_mac[lev])->array(mfi);
            const auto& wmac_fab = (m_w_mac[lev])->array(mfi);

            if(flags.getType(amrex::grow(bx, 0)) == FabType::covered)
            {
                m_u_mac[lev]->setVal(1.2345e300, ubx, 0, 1);
                m_v_mac[lev]->setVal(1.2345e300, vbx, 0, 1);
                m_w_mac[lev]->setVal(1.2345e300, wbx, 0, 1);
            }
            else if(flags.getType(amrex::grow(bx, 1)) == FabType::regular)
            {
                // No cut cells in tile + 1-cell witdh halo -> use non-eb routine
                AMREX_HOST_DEVICE_FOR_3D(ubx, i, j, k,
                {
                    // X-faces
                    Real upls     = ccvel_fab(i  ,j,k,0) - 0.5 * xslopes_fab(i  ,j,k,0);
                    Real umns     = ccvel_fab(i-1,j,k,0) + 0.5 * xslopes_fab(i-1,j,k,0);
                    if ( umns < 0.0 && upls > 0.0 )
                    {
                        umac_fab(i,j,k) = 0.0;
                    }
                    else
                    {
                        Real avg = 0.5 * ( upls + umns );
                        if (std::abs(avg) <  small_vel)
                            umac_fab(i,j,k) = 0.0;
                        else if (avg >= 0)
                            umac_fab(i,j,k) = umns;
                        else
                            umac_fab(i,j,k) = upls;
                    }
                });

                AMREX_HOST_DEVICE_FOR_3D(vbx, i, j, k,
                {
                    // Y-faces
                    Real upls     = ccvel_fab(i,j  ,k,1) - 0.5 * yslopes_fab(i,j  ,k,1);
                    Real umns     = ccvel_fab(i,j-1,k,1) + 0.5 * yslopes_fab(i,j-1,k,1);
                    if ( umns < 0.0 && upls > 0.0 )
                    {
                        vmac_fab(i,j,k) = 0.0;
                    }
                    else
                    {
                        Real avg = 0.5 * ( upls + umns );
                        if (std::abs(avg) <  small_vel)
                            vmac_fab(i,j,k) = 0.0;
                        else if (avg >= 0)
                            vmac_fab(i,j,k) = umns;
                        else
                            vmac_fab(i,j,k) = upls;
                    }
                });

                AMREX_HOST_DEVICE_FOR_3D(wbx, i, j, k,
                {
                    // Z-faces
                    Real upls     = ccvel_fab(i,j,k  ,2) - 0.5 * zslopes_fab(i,j,k  ,2);
                    Real umns     = ccvel_fab(i,j,k-1,2) + 0.5 * zslopes_fab(i,j,k-1,2);
                    if ( umns < 0.0 && upls > 0.0 )
                    {
                        wmac_fab(i,j,k) = 0.0;
                    }
                    else
                    {
                        Real avg = 0.5 * ( upls + umns );
                        if ( std::abs(avg) <  small_vel)
                            wmac_fab(i,j,k) = 0.0;
                        else if (avg >= 0)
                            wmac_fab(i,j,k) = umns;
                        else
                            wmac_fab(i,j,k) = upls;
                    }
                });

            }
            else
            {

                // Face-centered areas
                const auto& ax_fab = areafrac[0]->array(mfi);
                const auto& ay_fab = areafrac[1]->array(mfi);
                const auto& az_fab = areafrac[2]->array(mfi);

                // This FAB has cut cells
                AMREX_HOST_DEVICE_FOR_3D(ubx, i, j, k,
                {
                    // X-faces
                    if (ax_fab(i,j,k) > 0.0)
                    {
                        Real upls     = ccvel_fab(i  ,j,k,0) - 0.5 * xslopes_fab(i  ,j,k,0);
                        Real umns     = ccvel_fab(i-1,j,k,0) + 0.5 * xslopes_fab(i-1,j,k,0);
                        if ( umns < 0.0 && upls > 0.0 )
                        {
                            umac_fab(i,j,k) = 0.0;
                        }
                        else
                        {
                            Real avg = 0.5 * ( upls + umns );
                            if (std::abs(avg) <  small_vel)
                                umac_fab(i,j,k) = 0.0;
                            else if (avg >= 0)
                                umac_fab(i,j,k) = umns;
                            else
                                umac_fab(i,j,k) = upls;
                        }
                    }
                    else
                    {
                        umac_fab(i,j,k) = huge_vel;
                    }
                });

                AMREX_HOST_DEVICE_FOR_3D(vbx, i, j, k,
                {
                    // Y-faces
                    if (ay_fab(i,j,k) > 0.0)
                    {
                        Real upls     = ccvel_fab(i,j  ,k,1) - 0.5 * yslopes_fab(i,j  ,k,1);
                        Real umns     = ccvel_fab(i,j-1,k,1) + 0.5 * yslopes_fab(i,j-1,k,1);
                        if ( umns < 0.0 && upls > 0.0 )
                        {
                            vmac_fab(i,j,k) = 0.0;
                        }
                        else
                        {
                            Real avg = 0.5 * ( upls + umns );
                            if ( std::abs(avg) <  small_vel)
                                vmac_fab(i,j,k) = 0.0;
                            else if (avg >= 0)
                                vmac_fab(i,j,k) = umns;
                            else
                                vmac_fab(i,j,k) = upls;
                        }
                    }
                    else
                    {
                        vmac_fab(i,j,k) = huge_vel;
                    }
                });

                AMREX_HOST_DEVICE_FOR_3D(wbx, i, j, k,
                {
                    // Z-faces
                    if (az_fab(i,j,k) > 0.0)
                    {
                       Real upls     = ccvel_fab(i,j,k  ,2) - 0.5 * zslopes_fab(i,j,k  ,2);
                       Real umns     = ccvel_fab(i,j,k-1,2) + 0.5 * zslopes_fab(i,j,k-1,2);
                       if ( umns < 0.0 && upls > 0.0 )
                       {
                            wmac_fab(i,j,k) = 0.0;
                       }
                       else
                       {
                            Real avg = 0.5 * ( upls + umns );
                            if (std::abs(avg) <  small_vel)
                                wmac_fab(i,j,k) = 0.0;
                            else if (avg >= 0)
                                wmac_fab(i,j,k) = umns;
                            else
                                wmac_fab(i,j,k) = upls;
                       }
                    }
                    else
                    {
                        wmac_fab(i,j,k) = huge_vel;
                    }
                });

            } // Cut cells
        } // MFIter
    } // Levels
}
//
// Compute the slopes of each velocity component in all three directions.
//
void incflo::ComputeVelocitySlopes(int lev, MultiFab& Sborder)
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
			xslopes[lev]->setVal(1.2345e300, bx, 0, AMREX_SPACEDIM);
			yslopes[lev]->setVal(1.2345e300, bx, 0, AMREX_SPACEDIM);
			zslopes[lev]->setVal(1.2345e300, bx, 0, AMREX_SPACEDIM);
		}
		else
		{
            const auto&  vel_fab =      Sborder.array(mfi);
            const auto&   xs_fab = xslopes[lev]->array(mfi);
            const auto&   ys_fab = yslopes[lev]->array(mfi);
            const auto&   zs_fab = zslopes[lev]->array(mfi);

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
                   xs_fab(i,j,k,dir) = (du_xc       > 0.0) ? xslope : -xslope;

                   // Y direction
                   Real du_yl = 2.0*(vel_fab(i,j  ,k,dir) - vel_fab(i,j-1,k,dir));
                   Real du_yr = 2.0*(vel_fab(i,j+1,k,dir) - vel_fab(i,j  ,k,dir));
                   Real du_yc = 0.5*(vel_fab(i,j+1,k,dir) - vel_fab(i,j-1,k,dir));

                   Real yslope = amrex::min(std::abs(du_yl),std::abs(du_yc),std::abs(du_yr));
                   yslope            = (du_yr*du_yl > 0.0) ? yslope : 0.0;
                   ys_fab(i,j,k,dir) = (du_yc       > 0.0) ? yslope : -yslope;

                   // Z direction
                   Real du_zl = 2.0*(vel_fab(i,j,k  ,dir) - vel_fab(i,j,k-1,dir));
                   Real du_zr = 2.0*(vel_fab(i,j,k+1,dir) - vel_fab(i,j,k  ,dir));
                   Real du_zc = 0.5*(vel_fab(i,j,k+1,dir) - vel_fab(i,j,k-1,dir));

                   Real zslope = amrex::min(std::abs(du_zl),std::abs(du_zc),std::abs(du_zr));
                   zslope            = (du_zr*du_zl > 0.0) ? zslope : 0.0;
                   zs_fab(i,j,k,dir) = (du_zc       > 0.0) ? zslope : -zslope;
                });
			}
			else
			{
                const auto& flag_fab = flags.array();

                AMREX_HOST_DEVICE_FOR_4D(bx, ncomp, i, j, k, dir,
                {
                    if (flag_fab(i,j,k).isCovered())
                    {
                        xs_fab(i,j,k,dir) = 0.0;
                        ys_fab(i,j,k,dir) = 0.0;
                        zs_fab(i,j,k,dir) = 0.0;

                    } else {

                        // X direction
                        Real du_xl = (flag_fab(i-1,j,k).isCovered()) ? 0.0 :
                                     2.0*(vel_fab(i  ,j,k,dir) - vel_fab(i-1,j,k,dir));
                        Real du_xr = (flag_fab(i+1,j,k).isCovered()) ? 0.0 :
                                     2.0*(vel_fab(i+1,j,k,dir) - vel_fab(i  ,j,k,dir));
                        Real du_xc = 0.5*(vel_fab(i+1,j,k,dir) - vel_fab(i-1,j,k,dir));

                        Real xslope = amrex::min(std::abs(du_xl),std::abs(du_xc),std::abs(du_xr));
                        xslope            = (du_xr*du_xl > 0.0) ? xslope : 0.0;
                        xs_fab(i,j,k,dir) = (du_xc       > 0.0) ? xslope : -xslope;

                        // Y direction
                        Real du_yl = (flag_fab(i,j-1,k).isCovered()) ? 0.0 :
                                     2.0*(vel_fab(i,j  ,k,dir) - vel_fab(i,j-1,k,dir));
                        Real du_yr = (flag_fab(i,j+1,k).isCovered()) ? 0.0 :
                                     2.0*(vel_fab(i,j+1,k,dir) - vel_fab(i,j  ,k,dir));
                        Real du_yc = 0.5*(vel_fab(i,j+1,k,dir) - vel_fab(i,j-1,k,dir));

                        Real yslope = amrex::min(std::abs(du_yl),std::abs(du_yc),std::abs(du_yr));
                        yslope            = (du_yr*du_yl > 0.0) ? yslope : 0.0;
                        ys_fab(i,j,k,dir) = (du_yc       > 0.0) ? yslope : -yslope;

                        // Z direction
                        Real du_zl = (flag_fab(i,j,k-1).isCovered()) ? 0.0 :
                                     2.0*(vel_fab(i,j,k  ,dir) - vel_fab(i,j,k-1,dir));
                        Real du_zr = (flag_fab(i,j,k+1).isCovered()) ? 0.0 :
                                     2.0*(vel_fab(i,j,k+1,dir) - vel_fab(i,j,k  ,dir));
                        Real du_zc = 0.5*(vel_fab(i,j,k+1,dir) - vel_fab(i,j,k-1,dir));

                        Real zslope = amrex::min(std::abs(du_zl),std::abs(du_zc),std::abs(du_zr));
                        zslope          = (du_zr*du_zl > 0.0) ? zslope : 0.0;
                        zs_fab(i,j,k,dir) = (du_zc       > 0.0) ? zslope : -zslope;
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
                    xs_fab(i,j,k,dir) = (du_xc       > 0.0) ? xslope : -xslope;
                }
                if ( (i == domain.bigEnd(0)) && !flag_fab(i,j,k).isCovered() && ihi_ifab(i+1,j,k,0) == 20)
                {
                    Real du_xl = 2.0*(vel_fab(i  ,j,k,dir) - vel_fab(i-1,j,k,dir));
                    Real du_xr = 2.0*(vel_fab(i+1,j,k,dir) - vel_fab(i  ,j,k,dir));
                    Real du_xc = -(vel_fab(i-1,j,k,dir)+3.0*vel_fab(i,j,k,dir)-4.0*vel_fab(i+1,j,k,dir))/3.0;

                    Real xslope = amrex::min(std::abs(du_xl),std::abs(du_xc),std::abs(du_xr));
                    xslope            = (du_xr*du_xl > 0.0) ? xslope : 0.0;
                    xs_fab(i,j,k,dir) = (du_xc       > 0.0) ? xslope : -xslope;
                }

                if ( (j == domain.smallEnd(1)) && !flag_fab(i,j,k).isCovered() && jlo_ifab(i,j-1,k,0) == 20)
                {
                    Real du_yl = 2.0*(vel_fab(i,j  ,k,dir) - vel_fab(i,j-1,k,dir));
                    Real du_yr = 2.0*(vel_fab(i,j+1,k,dir) - vel_fab(i,j  ,k,dir));
                    Real du_yc = (vel_fab(i,j+1,k,dir)+3.0*vel_fab(i,j,k,dir)-4.0*vel_fab(i,j-1,k,dir))/3.0;

                    Real yslope = amrex::min(std::abs(du_yl),std::abs(du_yc),std::abs(du_yr));
                    yslope            = (du_yr*du_yl > 0.0) ? yslope : 0.0;
                    ys_fab(i,j,k,dir) = (du_yc       > 0.0) ? yslope : -yslope;
                }
                if ( (j == domain.bigEnd(1)) && !flag_fab(i,j,k).isCovered() && jhi_ifab(i,j+1,k,0) == 20)
                {
                    Real du_yl = 2.0*(vel_fab(i,j  ,k,dir) - vel_fab(i,j-1,k,dir));
                    Real du_yr = 2.0*(vel_fab(i,j+1,k,dir) - vel_fab(i,j  ,k,dir));
                    Real du_yc = -(vel_fab(i,j-1,k,dir)+3.0*vel_fab(i,j,k,dir)-4.0*vel_fab(i,j+1,k,dir))/3.0;

                    Real yslope = amrex::min(std::abs(du_yl),std::abs(du_yc),std::abs(du_yr));
                    yslope            = (du_yr*du_yl > 0.0) ? yslope : 0.0;
                    ys_fab(i,j,k,dir) = (du_yc       > 0.0) ? yslope : -yslope;
                }

                if ( (k == domain.smallEnd(2)) && !flag_fab(i,j,k).isCovered() && klo_ifab(i,j,k-1,0) == 20)
                {
                    Real du_zl = 2.0*(vel_fab(i,j,k  ,dir) - vel_fab(i,j,k-1,dir));
                    Real du_zr = 2.0*(vel_fab(i,j,k+1,dir) - vel_fab(i,j,k  ,dir));
                    Real du_zc = (vel_fab(i,j,k+1,dir)+3.0*vel_fab(i,j,k,dir)-4.0*vel_fab(i,j,k-1,dir))/3.0;

                    Real zslope = amrex::min(std::abs(du_zl),std::abs(du_zc),std::abs(du_zr));
                    zslope            = (du_zr*du_zl > 0.0) ? zslope : 0.0;
                    zs_fab(i,j,k,dir) = (du_zc       > 0.0) ? zslope : -zslope;
                }
                if ( (k == domain.bigEnd(2)) && !flag_fab(i,j,k).isCovered() && khi_ifab(i,j,k+1,0) == 20)
                {
                    Real du_zl = 2.0*(vel_fab(i,j,k  ,dir) - vel_fab(i,j,k-1,dir));
                    Real du_zr = 2.0*(vel_fab(i,j,k+1,dir) - vel_fab(i,j,k  ,dir));
                    Real du_zc = -(vel_fab(i,j,k-1,dir)+3.0*vel_fab(i,j,k,dir)-4.0*vel_fab(i,j,k+1,dir))/3.0;

                    Real zslope = amrex::min(std::abs(du_zl),std::abs(du_zc),std::abs(du_zr));
                    zslope            = (du_zr*du_zl > 0.0) ? zslope : 0.0;
                    zs_fab(i,j,k,dir) = (du_zc       > 0.0) ? zslope : -zslope;
                }
            });
		}
	}

	xslopes[lev]->FillBoundary(geom[lev].periodicity());
	yslopes[lev]->FillBoundary(geom[lev].periodicity());
	zslopes[lev]->FillBoundary(geom[lev].periodicity());
}

//
// Compute the three components of the convection term
//
void
incflo::incflo_compute_ugradu(Box& bx,
                              Vector< std::unique_ptr<MultiFab> >& conv, 
                              Vector< std::unique_ptr<MultiFab> >& vel,
                              Vector< std::unique_ptr<MultiFab> >& u_mac,
                              Vector< std::unique_ptr<MultiFab> >& v_mac,
                              Vector< std::unique_ptr<MultiFab> >& w_mac,
                              MFIter* mfi,
                              Box& domain,
                              const int lev)
{
  const Real* dx = geom[lev].CellSize();
  const amrex::Dim3 dom_low = amrex::lbound(domain);
  const amrex::Dim3 dom_high = amrex::ubound(domain);

  Array4<Real> const& ugradu = conv[lev]->array(*mfi); 
  
  Array4<Real> const& velocity = vel[lev]->array(*mfi);
  
  Array4<Real> const& u = u_mac[lev]->array(*mfi);
  Array4<Real> const& v = v_mac[lev]->array(*mfi);
  Array4<Real> const& w = w_mac[lev]->array(*mfi);

  Array4<Real> const& x_slopes = xslopes[lev]->array(*mfi);
  Array4<Real> const& y_slopes = yslopes[lev]->array(*mfi);
  Array4<Real> const& z_slopes = zslopes[lev]->array(*mfi);

  Array4<int> const& bc_ilo_type = bc_ilo[lev]->array();
  Array4<int> const& bc_ihi_type = bc_ihi[lev]->array();
  Array4<int> const& bc_jlo_type = bc_jlo[lev]->array();
  Array4<int> const& bc_jhi_type = bc_jhi[lev]->array();
  Array4<int> const& bc_klo_type = bc_klo[lev]->array();
  Array4<int> const& bc_khi_type = bc_khi[lev]->array();

  const Real i_dx(1/dx[0]), i_dy(1/dx[1]), i_dz(1/dx[2]);

  // Vectorize the boundary conditions list in order to use it in lambda
  // functions
  const GpuArray<int, 3> bc_types =
    {bc_list.get_minf(), bc_list.get_pinf(), bc_list.get_pout()};

  AMREX_HOST_DEVICE_FOR_3D(bx, i, j, k,
  {
    Real u_w(0); Real v_w(0); Real w_w(0);
    Real u_e(0); Real v_e(0); Real w_e(0);
    Real u_s(0); Real v_s(0); Real w_s(0);
    Real u_n(0); Real v_n(0); Real w_n(0);
    Real u_b(0); Real v_b(0); Real w_b(0);
    Real u_t(0); Real v_t(0); Real w_t(0);
    Real umns(0); Real vmns(0); Real wmns(0);
    Real upls(0); Real vpls(0); Real wpls(0);

    //
    // West face
    //
    // In the case of MINF       we are using the prescribed Dirichlet value
    // In the case of PINF, POUT we are using the upwind value
    if((i == dom_low.x) and
     ugradu_aux::is_equal_to_any(bc_ilo_type(dom_low.x-1,j,k,0),
                                 bc_types.data(), bc_types.size()))
    {
      u_w = velocity(i-1,j,k,0);
      v_w = velocity(i-1,j,k,1);
      w_w = velocity(i-1,j,k,2);
    }
    else {
      upls = velocity(i  ,j,k,0) - .5*x_slopes(i  ,j,k,0);
      umns = velocity(i-1,j,k,0) + .5*x_slopes(i-1,j,k,0);
      vpls = velocity(i  ,j,k,1) - .5*x_slopes(i  ,j,k,1);
      vmns = velocity(i-1,j,k,1) + .5*x_slopes(i-1,j,k,1);
      wpls = velocity(i  ,j,k,2) - .5*x_slopes(i  ,j,k,2);
      wmns = velocity(i-1,j,k,2) + .5*x_slopes(i-1,j,k,2);

      u_w = upwind( umns, upls, u(i,j,k) );
      v_w = upwind( vmns, vpls, u(i,j,k) );
      w_w = upwind( wmns, wpls, u(i,j,k) );
    }

    //
    // East face
    //
    // In the case of MINF       we are using the prescribed Dirichlet value
    // In the case of PINF, POUT we are using the upwind value
    if((i == dom_high.x) and
     ugradu_aux::is_equal_to_any(bc_ihi_type(dom_high.x+1,j,k,0),
                                 bc_types.data(), bc_types.size()))
    {
      u_e = velocity(i+1,j,k,0);
      v_e = velocity(i+1,j,k,1);
      w_e = velocity(i+1,j,k,2);
    }
    else {
      upls = velocity(i+1,j,k,0) - .5*x_slopes(i+1,j,k,0);
      umns = velocity(i  ,j,k,0) + .5*x_slopes(i  ,j,k,0);
      vpls = velocity(i+1,j,k,1) - .5*x_slopes(i+1,j,k,1);
      vmns = velocity(i  ,j,k,1) + .5*x_slopes(i  ,j,k,1);
      wpls = velocity(i+1,j,k,2) - .5*x_slopes(i+1,j,k,2);
      wmns = velocity(i  ,j,k,2) + .5*x_slopes(i  ,j,k,2);

      u_e = upwind( umns, upls, u(i+1,j,k) );
      v_e = upwind( vmns, vpls, u(i+1,j,k) );
      w_e = upwind( wmns, wpls, u(i+1,j,k) );
    }

    //
    // South face
    //
    // In the case of MINF       we are using the prescribed Dirichlet value
    // In the case of PINF, POUT we are using the upwind value
    if((j == dom_low.y) and
     ugradu_aux::is_equal_to_any(bc_jlo_type(i,dom_low.y-1,k,0),
                                 bc_types.data(), bc_types.size()))
    {
      u_s = velocity(i,j-1,k,0);
      v_s = velocity(i,j-1,k,1);
      w_s = velocity(i,j-1,k,2);
    }
    else {
      upls = velocity(i,j  ,k,0) - .5*y_slopes(i,j  ,k,0);
      umns = velocity(i,j-1,k,0) + .5*y_slopes(i,j-1,k,0);
      vpls = velocity(i,j  ,k,1) - .5*y_slopes(i,j  ,k,1);
      vmns = velocity(i,j-1,k,1) + .5*y_slopes(i,j-1,k,1);
      wpls = velocity(i,j  ,k,2) - .5*y_slopes(i,j  ,k,2);
      wmns = velocity(i,j-1,k,2) + .5*y_slopes(i,j-1,k,2);

      u_s = upwind( umns, upls, v(i,j,k) );
      v_s = upwind( vmns, vpls, v(i,j,k) );
      w_s = upwind( wmns, wpls, v(i,j,k) );
    }

    //
    // North face
    //
    // In the case of MINF       we are using the prescribed Dirichlet value
    // In the case of PINF, POUT we are using the upwind value
    if((j == dom_high.y) and
     ugradu_aux::is_equal_to_any(bc_jhi_type(i,dom_high.y+1,k,0),
                                 bc_types.data(), bc_types.size()))
    {
      u_n = velocity(i,j+1,k,0);
      v_n = velocity(i,j+1,k,1);
      w_n = velocity(i,j+1,k,2);
    }
    else {
      upls = velocity(i,j+1,k,0) - .5*y_slopes(i,j+1,k,0);
      umns = velocity(i,j  ,k,0) + .5*y_slopes(i,j  ,k,0);
      vpls = velocity(i,j+1,k,1) - .5*y_slopes(i,j+1,k,1);
      vmns = velocity(i,j  ,k,1) + .5*y_slopes(i,j  ,k,1);
      wpls = velocity(i,j+1,k,2) - .5*y_slopes(i,j+1,k,2);
      wmns = velocity(i,j  ,k,2) + .5*y_slopes(i,j  ,k,2);

      u_n = upwind( umns, upls, v(i,j+1,k) );
      v_n = upwind( vmns, vpls, v(i,j+1,k) );
      w_n = upwind( wmns, wpls, v(i,j+1,k) );
    }

    //
    // Bottom face
    //
    // In the case of MINF       we are using the prescribed Dirichlet value
    // In the case of PINF, POUT we are using the upwind value
    if((k == dom_low.z) and
     ugradu_aux::is_equal_to_any(bc_klo_type(i,j,dom_low.z-1,0),
                                 bc_types.data(), bc_types.size()))
    {
      u_b = velocity(i,j,k-1,0);
      v_b = velocity(i,j,k-1,1);
      w_b = velocity(i,j,k-1,2);
    }
    else {
      upls = velocity(i,j,k  ,0) - .5*z_slopes(i,j,k  ,0);
      umns = velocity(i,j,k-1,0) + .5*z_slopes(i,j,k-1,0);
      vpls = velocity(i,j,k  ,1) - .5*z_slopes(i,j,k  ,1);
      vmns = velocity(i,j,k-1,1) + .5*z_slopes(i,j,k-1,1);
      wpls = velocity(i,j,k  ,2) - .5*z_slopes(i,j,k  ,2);
      wmns = velocity(i,j,k-1,2) + .5*z_slopes(i,j,k-1,2);

      u_b = upwind( umns, upls, w(i,j,k) );
      v_b = upwind( vmns, vpls, w(i,j,k) );
      w_b = upwind( wmns, wpls, w(i,j,k) );
    }

    //
    // Top face
    //
    // In the case of MINF       we are using the prescribed Dirichlet value
    // In the case of PINF, POUT we are using the upwind value
    if((k == dom_high.z) and
     ugradu_aux::is_equal_to_any(bc_khi_type(i,j,dom_high.z+1,0),
                                 bc_types.data(), bc_types.size()))
    {
      u_t = velocity(i,j,k+1,0);
      v_t = velocity(i,j,k+1,1);
      w_t = velocity(i,j,k+1,2);
    }
    else {
      upls = velocity(i,j,k+1,0) - .5*z_slopes(i,j,k+1,0);
      umns = velocity(i,j,k  ,0) + .5*z_slopes(i,j,k  ,0);
      vpls = velocity(i,j,k+1,1) - .5*z_slopes(i,j,k+1,1);
      vmns = velocity(i,j,k  ,1) + .5*z_slopes(i,j,k  ,1);
      wpls = velocity(i,j,k+1,2) - .5*z_slopes(i,j,k+1,2);
      wmns = velocity(i,j,k  ,2) + .5*z_slopes(i,j,k  ,2);

      u_t = upwind( umns, upls, w(i,j,k+1) );
      v_t = upwind( vmns, vpls, w(i,j,k+1) );
      w_t = upwind( wmns, wpls, w(i,j,k+1) );
    }

    // Define the convective terms -- conservatively
    //   ugradu = ( div(u^MAC u^cc) - u^cc div(u^MAC) ) 
    Real u_hi_x = u(i+1,j,k);
    Real u_lo_x = u(i  ,j,k);
    Real v_hi_y = v(i,j+1,k);
    Real v_lo_y = v(i,j  ,k);
    Real w_hi_z = w(i,j,k+1);
    Real w_lo_z = w(i,j,k  );

    Real divumac = (u_hi_x - u_lo_x) * i_dx +
                   (v_hi_y - v_lo_y) * i_dy + 
                   (w_hi_z - w_lo_z) * i_dz;

    ugradu(i,j,k,0) = (u_hi_x*u_e - u_lo_x*u_w) * i_dx +
                      (v_hi_y*u_n - v_lo_y*u_s) * i_dy +
                      (w_hi_z*u_t - w_lo_z*u_b) * i_dz -
                      velocity(i,j,k,0)*divumac;

    ugradu(i,j,k,1) = (u_hi_x*v_e - u_lo_x*v_w) * i_dx +
                      (v_hi_y*v_n - v_lo_y*v_s) * i_dy +
                      (w_hi_z*v_t - w_lo_z*v_b) * i_dz -
                      velocity(i,j,k,1)*divumac;

    ugradu(i,j,k,2) = (u_hi_x*w_e - u_lo_x*w_w) * i_dx +
                      (v_hi_y*w_n - v_lo_y*w_s) * i_dy +
                      (w_hi_z*w_t - w_lo_z*w_b) * i_dz -
                      velocity(i,j,k,2)*divumac;

    //
    // Return the negative
    //
    const Real coefficient(-1.);
    ugradu(i,j,k,0) *= coefficient; 
    ugradu(i,j,k,1) *= coefficient; 
    ugradu(i,j,k,2) *= coefficient; 
  });

#ifdef AMREX_USE_CUDA
  Gpu::Device::synchronize();
#endif
}
//
// Compute the three components of the convection term when we have embedded
// boundaries
//
void
incflo::incflo_compute_ugradu_eb(Box& bx,
                                 Vector< std::unique_ptr<MultiFab> >& conv, 
                                 Vector< std::unique_ptr<MultiFab> >& vel,
                                 Vector< std::unique_ptr<MultiFab> >& u_mac,
                                 Vector< std::unique_ptr<MultiFab> >& v_mac,
                                 Vector< std::unique_ptr<MultiFab> >& w_mac,
                                 MFIter* mfi,
                                 Array<const MultiCutFab*,AMREX_SPACEDIM>& areafrac,
                                 Array<const MultiCutFab*,AMREX_SPACEDIM>& facecent,
                                 const MultiFab* volfrac,
                                 const MultiCutFab* bndrycent,
                                 Box& domain,
                                 const EBCellFlagFab& flags,
                                 const int lev)
{
  AMREX_ASSERT_WITH_MESSAGE(nghost >= 4, "Compute divop_conv(): ng must be >= 4");

  const Real* dx = geom[lev].CellSize();
  const amrex::Dim3 dom_low = amrex::lbound(domain);
  const amrex::Dim3 dom_high = amrex::ubound(domain);

  Array4<Real> const& ugradu   = conv[lev]->array(*mfi);
  Array4<Real> const& velocity =  vel[lev]->array(*mfi);

  Array4<const Real> const& areafrac_x = areafrac[0]->array(*mfi);
  Array4<const Real> const& areafrac_y = areafrac[1]->array(*mfi);
  Array4<const Real> const& areafrac_z = areafrac[2]->array(*mfi);

  Array4<Real> const& u = u_mac[lev]->array(*mfi);
  Array4<Real> const& v = v_mac[lev]->array(*mfi);
  Array4<Real> const& w = w_mac[lev]->array(*mfi);

  Array4<Real> const& x_slopes = xslopes[lev]->array(*mfi);
  Array4<Real> const& y_slopes = yslopes[lev]->array(*mfi);
  Array4<Real> const& z_slopes = zslopes[lev]->array(*mfi);

  Array4<int> const& bc_ilo_type = bc_ilo[lev]->array();
  Array4<int> const& bc_ihi_type = bc_ihi[lev]->array();
  Array4<int> const& bc_jlo_type = bc_jlo[lev]->array();
  Array4<int> const& bc_jhi_type = bc_jhi[lev]->array();
  Array4<int> const& bc_klo_type = bc_klo[lev]->array();
  Array4<int> const& bc_khi_type = bc_khi[lev]->array();

  // Number of Halo layers
  const int nh(3);

  const Box ubx = amrex::surroundingNodes(amrex::grow(bx,nh),0);
  const Box vbx = amrex::surroundingNodes(amrex::grow(bx,nh),1);
  const Box wbx = amrex::surroundingNodes(amrex::grow(bx,nh),2);

  const int ncomp(3);
  
  FArrayBox fxfab(ubx, ncomp);
  FArrayBox fyfab(vbx, ncomp);
  FArrayBox fzfab(wbx, ncomp);

  Array4<Real> const& fx = fxfab.array();
  Array4<Real> const& fy = fyfab.array();
  Array4<Real> const& fz = fzfab.array();

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
  AMREX_HOST_DEVICE_FOR_4D(ubx, ncomp, i, j, k, n,
  {
    Real u_face(0);
    Real upls(0); Real umns(0);

    if( areafrac_x(i,j,k) > 0 ) {
      if( i <= dom_low.x and
       ugradu_aux::is_equal_to_any(bc_ilo_type(dom_low.x-1,j,k,0),
                                   bc_types.data(), bc_types.size()))
      {
        u_face = velocity(dom_low.x-1,j,k,n);
      }
      else if( i >= dom_high.x+1 and
       ugradu_aux::is_equal_to_any(bc_ihi_type(dom_high.x+1,j,k,0),
                                   bc_types.data(), bc_types.size()))
      {
        u_face = velocity(dom_high.x+1,j,k,n);
      }
      else {
        upls = velocity(i  ,j,k,n) - .5*x_slopes(i  ,j,k,n);
        umns = velocity(i-1,j,k,n) + .5*x_slopes(i-1,j,k,n);

        u_face = upwind( umns, upls, u(i,j,k) );
      }
    }
    else {
      u_face = my_huge; 
    }
    fx(i,j,k,n) = u(i,j,k) * u_face;
  });

  //
  // ===================== Y =====================
  //
  AMREX_HOST_DEVICE_FOR_4D(vbx, ncomp, i, j, k, n,
  {
    Real v_face(0);
    Real vpls(0); Real vmns(0);

    if( areafrac_y(i,j,k) > 0 ) {
      if( j <= dom_low.y and
       ugradu_aux::is_equal_to_any(bc_jlo_type(i,dom_low.y-1,k,0),
                                   bc_types.data(), bc_types.size()))
      {
        v_face = velocity(i,dom_low.y-1,k,n);
      }
      else if( j >= dom_high.y+1 and
       ugradu_aux::is_equal_to_any(bc_jhi_type(i,dom_high.y+1,k,0),
                                   bc_types.data(), bc_types.size()))
      {
        v_face = velocity(i,dom_high.y+1,k,n);
      }
      else {
        vpls = velocity(i,j  ,k,n) - .5*y_slopes(i,j  ,k,n);
        vmns = velocity(i,j-1,k,n) + .5*y_slopes(i,j-1,k,n);

        v_face = upwind( vmns, vpls, v(i,j,k) );
      }
    }
    else {
      v_face = my_huge;
    }
    fy(i,j,k,n) = v(i,j,k) * v_face;
  });

  //
  // ===================== Z =====================
  //
  AMREX_HOST_DEVICE_FOR_4D(wbx, ncomp, i, j, k, n,
  {
    Real w_face(0);
    Real wpls(0); Real wmns(0);

    if( areafrac_z(i,j,k) > 0 ) {
      if( k <= dom_low.z and
       ugradu_aux::is_equal_to_any(bc_klo_type(i,j,dom_low.z-1,0),
                                   bc_types.data(), bc_types.size()))
      {
        w_face = velocity(i,j,dom_low.z-1,n);
      }
      else if( k >= dom_high.z+1 and
       ugradu_aux::is_equal_to_any(bc_khi_type(i,j,dom_high.z+1,0),
                                   bc_types.data(), bc_types.size()))
      {
        w_face = velocity(i,j,dom_high.z+1,n);
      }
      else {
        wpls = velocity(i,j,k  ,n) - .5*z_slopes(i,j,k  ,n);
        wmns = velocity(i,j,k-1,n) + .5*z_slopes(i,j,k-1,n);

        w_face = upwind( wmns, wpls, w(i,j,k) );
      }
    }
    else {
      w_face = my_huge;
    }
    fz(i,j,k,n) = w(i,j,k) * w_face;
  });

#ifdef AMREX_USE_CUDA
  Gpu::Device::synchronize();
#endif

  const int cyclic_x = geom[0].isPeriodic(0) ? 1 : 0;
  const int cyclic_y = geom[0].isPeriodic(1) ? 1 : 0;
  const int cyclic_z = geom[0].isPeriodic(2) ? 1 : 0;

  // Compute div(tau) with EB algorithm
  compute_divop_conv(bx, *conv[lev], mfi, fxfab, fyfab, fzfab, 
                     areafrac, facecent, flags, volfrac, bndrycent, domain,
                     cyclic_x, cyclic_y, cyclic_z, dx);

  AMREX_HOST_DEVICE_FOR_3D(bx, i, j, k,
  {
    const Real coefficient(-1.);
    ugradu(i,j,k,0) *= coefficient; 
    ugradu(i,j,k,1) *= coefficient; 
    ugradu(i,j,k,2) *= coefficient; 
  });

#ifdef AMREX_USE_CUDA
  Gpu::Device::synchronize();
#endif
}
