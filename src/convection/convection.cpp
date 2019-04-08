#include <AMReX_EBMultiFabUtil.H>

#include <incflo.H>
#include <mac_F.H>
#include <convection_F.H>

//
// Compute acc using the vel passed in
//
void incflo::ComputeUGradU(Vector<std::unique_ptr<MultiFab>>& conv_in,
                           Vector<std::unique_ptr<MultiFab>>& vel_in, 
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
                conv_in[lev]->setVal(1.2345e300, bx, 0, 3);
            }
            else
            {
                // No cut cells in tile + nghost-cell witdh halo -> use non-eb routine
                if(flags.getType(amrex::grow(bx, nghost)) == FabType::regular)
                {
                    compute_ugradu(BL_TO_FORTRAN_BOX(bx),
                                   BL_TO_FORTRAN_ANYD((*conv_in[lev])[mfi]),
                                   BL_TO_FORTRAN_ANYD((*vel_in[lev])[mfi]),
                                   BL_TO_FORTRAN_ANYD((*m_u_mac[lev])[mfi]),
                                   BL_TO_FORTRAN_ANYD((*m_v_mac[lev])[mfi]),
                                   BL_TO_FORTRAN_ANYD((*m_w_mac[lev])[mfi]),
                                   (*xslopes[lev])[mfi].dataPtr(),
                                   (*yslopes[lev])[mfi].dataPtr(),
                                   BL_TO_FORTRAN_ANYD((*zslopes[lev])[mfi]),
                                   domain.loVect(),
                                   domain.hiVect(),
                                   bc_ilo[lev]->dataPtr(),
                                   bc_ihi[lev]->dataPtr(),
                                   bc_jlo[lev]->dataPtr(),
                                   bc_jhi[lev]->dataPtr(),
                                   bc_klo[lev]->dataPtr(),
                                   bc_khi[lev]->dataPtr(),
                                   geom[lev].CellSize(),
                                   &nghost);
                }
                else
                {
                    compute_ugradu_eb(BL_TO_FORTRAN_BOX(bx),
                                      BL_TO_FORTRAN_ANYD((*conv_in[lev])[mfi]),
                                      BL_TO_FORTRAN_ANYD((*vel_in[lev])[mfi]),
                                      BL_TO_FORTRAN_ANYD((*m_u_mac[lev])[mfi]),
                                      BL_TO_FORTRAN_ANYD((*m_v_mac[lev])[mfi]),
                                      BL_TO_FORTRAN_ANYD((*m_w_mac[lev])[mfi]),
                                      BL_TO_FORTRAN_ANYD((*areafrac[0])[mfi]),
                                      BL_TO_FORTRAN_ANYD((*areafrac[1])[mfi]),
                                      BL_TO_FORTRAN_ANYD((*areafrac[2])[mfi]),
                                      BL_TO_FORTRAN_ANYD((*facecent[0])[mfi]),
                                      BL_TO_FORTRAN_ANYD((*facecent[1])[mfi]),
                                      BL_TO_FORTRAN_ANYD((*facecent[2])[mfi]),
                                      BL_TO_FORTRAN_ANYD(flags),
                                      BL_TO_FORTRAN_ANYD((*volfrac)[mfi]),
                                      BL_TO_FORTRAN_ANYD((*bndrycent)[mfi]),
                                      (*xslopes[lev])[mfi].dataPtr(),
                                      (*yslopes[lev])[mfi].dataPtr(),
                                      BL_TO_FORTRAN_ANYD((*zslopes[lev])[mfi]),
                                      domain.loVect(),
                                      domain.hiVect(),
                                      bc_ilo[lev]->dataPtr(),
                                      bc_ihi[lev]->dataPtr(),
                                      bc_jlo[lev]->dataPtr(),
                                      bc_jhi[lev]->dataPtr(),
                                      bc_klo[lev]->dataPtr(),
                                      bc_khi[lev]->dataPtr(),
                                      geom[lev].CellSize(),
                                      &nghost);
                }
            }
        }
	}
}

//
// Compute the slopes of each velocity component in the
// three directions.
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
			xslopes[lev]->setVal(1.2345e300, bx, 0, 3);
			yslopes[lev]->setVal(1.2345e300, bx, 0, 3);
			zslopes[lev]->setVal(1.2345e300, bx, 0, 3);
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
                amrex::ParallelFor(bx, ncomp,
                      [=] (int i, int j, int k, int dir)
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

                amrex::ParallelFor(bx, ncomp,
                      [=] (int i, int j, int k, int dir)
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

            correct_slopes(BL_TO_FORTRAN_BOX(bx),
                           BL_TO_FORTRAN_ANYD(Sborder[mfi]),
                           (*xslopes[lev])[mfi].dataPtr(),
                           (*yslopes[lev])[mfi].dataPtr(),
                           BL_TO_FORTRAN_ANYD((*zslopes[lev])[mfi]),
                           BL_TO_FORTRAN_ANYD(flags),
                           domain.loVect(), domain.hiVect(),
                           bc_ilo[lev]->dataPtr(), bc_ihi[lev]->dataPtr(),
                           bc_jlo[lev]->dataPtr(), bc_jhi[lev]->dataPtr(),
                           bc_klo[lev]->dataPtr(), bc_khi[lev]->dataPtr(),
                           &nghost);
		}
	}

	xslopes[lev]->FillBoundary(geom[lev].periodicity());
	yslopes[lev]->FillBoundary(geom[lev].periodicity());
	zslopes[lev]->FillBoundary(geom[lev].periodicity());
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

            if(flags.getType(amrex::grow(bx, 0)) == FabType::covered)
            {
                m_u_mac[lev]->setVal(1.2345e300, ubx, 0, 1);
                m_v_mac[lev]->setVal(1.2345e300, vbx, 0, 1);
                m_w_mac[lev]->setVal(1.2345e300, wbx, 0, 1);
            }
            else
            {
                // No cut cells in tile + 1-cell witdh halo -> use non-eb routine
                if(flags.getType(amrex::grow(bx, 1)) == FabType::regular)
                {
                    compute_velocity_at_faces(BL_TO_FORTRAN_BOX(bx),
                                              BL_TO_FORTRAN_ANYD((*m_u_mac[lev])[mfi]),
                                              BL_TO_FORTRAN_ANYD((*m_v_mac[lev])[mfi]),
                                              BL_TO_FORTRAN_ANYD((*m_w_mac[lev])[mfi]),
                                              BL_TO_FORTRAN_ANYD((*vel_in[lev])[mfi]),
                                              BL_TO_FORTRAN_ANYD((*xslopes[lev])[mfi]),
                                              (*yslopes[lev])[mfi].dataPtr(),
                                              (*zslopes[lev])[mfi].dataPtr(),
                                              bc_ilo[lev]->dataPtr(),
                                              bc_ihi[lev]->dataPtr(),
                                              bc_jlo[lev]->dataPtr(),
                                              bc_jhi[lev]->dataPtr(),
                                              bc_klo[lev]->dataPtr(),
                                              bc_khi[lev]->dataPtr(),
                                              &nghost,
                                              domain.loVect(),
                                              domain.hiVect());
                }
                else
                {
                    compute_velocity_at_x_faces_eb(BL_TO_FORTRAN_BOX(ubx),
                                                   BL_TO_FORTRAN_ANYD((*m_u_mac[lev])[mfi]),
                                                   BL_TO_FORTRAN_ANYD((*vel_in[lev])[mfi]),
                                                   BL_TO_FORTRAN_ANYD((*xslopes[lev])[mfi]),
                                                   BL_TO_FORTRAN_ANYD((*areafrac[0])[mfi]),
                                                   BL_TO_FORTRAN_ANYD(flags),
                                                   &nghost,
                                                   domain.loVect(),
                                                   domain.hiVect());

                    compute_velocity_at_y_faces_eb(BL_TO_FORTRAN_BOX(vbx),
                                                   BL_TO_FORTRAN_ANYD((*m_v_mac[lev])[mfi]),
                                                   BL_TO_FORTRAN_ANYD((*vel_in[lev])[mfi]),
                                                   BL_TO_FORTRAN_ANYD((*yslopes[lev])[mfi]),
                                                   BL_TO_FORTRAN_ANYD((*areafrac[1])[mfi]),
                                                   BL_TO_FORTRAN_ANYD(flags),
                                                   &nghost,
                                                   domain.loVect(),
                                                   domain.hiVect());

                    compute_velocity_at_z_faces_eb(BL_TO_FORTRAN_BOX(wbx),
                                                   BL_TO_FORTRAN_ANYD((*m_w_mac[lev])[mfi]),
                                                   BL_TO_FORTRAN_ANYD((*vel_in[lev])[mfi]),
                                                   BL_TO_FORTRAN_ANYD((*zslopes[lev])[mfi]),
                                                   BL_TO_FORTRAN_ANYD((*areafrac[2])[mfi]),
                                                   BL_TO_FORTRAN_ANYD(flags),
                                                   &nghost,
                                                   domain.loVect(),
                                                   domain.hiVect());
                }
            }
        }
    }
}
