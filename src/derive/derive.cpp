#include <AMReX_Box.H>

#include <incflo.H>
#include <derive_F.H>
#include <incflo_proj_F.H>

void incflo::ComputeDivU(Real time_in)
{
    int extrap_dir_bcs = 0;
    incflo_set_velocity_bcs(time_in, vel, extrap_dir_bcs);

    nodal_projector->computeRHS(GetVecOfPtrs(divu),GetVecOfPtrs(vel));
}

void incflo::ComputeStrainrate(Real time_in)
{
    BL_PROFILE("incflo::ComputeStrainrate");

    for(int lev = 0; lev <= finest_level; lev++)
    {
        Box domain(geom[lev].Domain());
        Real idx = 1.0 / geom[lev].CellSize()[0];
        Real idy = 1.0 / geom[lev].CellSize()[1];
        Real idz = 1.0 / geom[lev].CellSize()[2];

        // State with ghost cells
        MultiFab Sborder(grids[lev], dmap[lev], vel[lev]->nComp(), nghost,
                         MFInfo(), *ebfactory[lev]);
        FillPatchVel(lev, time_in, Sborder);

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(Sborder, TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            // Tilebox
            Box bx = mfi.tilebox();

            // This is to check efficiently if this tile contains any eb stuff
            const EBFArrayBox& vel_fab = static_cast<EBFArrayBox const&>(Sborder[mfi]);
            const EBCellFlagFab& flags = vel_fab.getEBCellFlagFab();

            // Cell-centered velocity
            const auto& ccvel_fab = Sborder.array(mfi);

            // Cell-centered strain-rate magnitude
            const auto& sr_fab = strainrate[lev]->array(mfi);

            if (flags.getType(amrex::grow(bx, 0)) == FabType::covered)
            {
                (*strainrate[lev])[mfi].setVal(1.2345e200, bx);
            }
            else if(flags.getType(amrex::grow(bx, 1)) == FabType::regular)
            {
                // No cut cells in tile + 1-cell witdh halo -> use non-eb routine
                AMREX_FOR_3D(bx, i, j, k,
                {
                    Real ux = 0.5 * (ccvel_fab(i+1,j,k,0) - ccvel_fab(i-1,j,k,0)) * idx;
                    Real vx = 0.5 * (ccvel_fab(i+1,j,k,1) - ccvel_fab(i-1,j,k,1)) * idx;
                    Real wx = 0.5 * (ccvel_fab(i+1,j,k,2) - ccvel_fab(i-1,j,k,2)) * idx;

                    Real uy = 0.5 * (ccvel_fab(i,j+1,k,0) - ccvel_fab(i,j-1,k,0)) * idy;
                    Real vy = 0.5 * (ccvel_fab(i,j+1,k,1) - ccvel_fab(i,j-1,k,1)) * idy;
                    Real wy = 0.5 * (ccvel_fab(i,j+1,k,2) - ccvel_fab(i,j-1,k,2)) * idy;

                    Real uz = 0.5 * (ccvel_fab(i,j,k+1,0) - ccvel_fab(i,j,k-1,0)) * idz;
                    Real vz = 0.5 * (ccvel_fab(i,j,k+1,1) - ccvel_fab(i,j,k-1,1)) * idz;
                    Real wz = 0.5 * (ccvel_fab(i,j,k+1,2) - ccvel_fab(i,j,k-1,2)) * idz;

                    // Include the factor half here rather than in each of the above
                    sr_fab(i,j,k) = sqrt(2.0 * pow(ux, 2) + 2.0 * pow(vy, 2) + 2.0 * pow(wz, 2)
                            + pow(uy + vx, 2) + pow(vz + wy, 2) + pow(wx + uz, 2));
                });
            }
            else
            {
                // Cut cells present -> use EB routine!
                const auto& flag_fab = flags.array();
                Real c0 = -1.5;
                Real c1 = 2.0;
                Real c2 = -0.5;

                AMREX_FOR_3D(bx, i, j, k,
                {
                    if (flag_fab(i,j,k).isCovered())
                    {
                        // Don't compute strainrate in cut cells
                        sr_fab(i,j,k) = 1.2345e200;
                    }
                    else
                    {
                        Real ux(0); Real vx(0); Real wx(0);
                        Real uy(0); Real vy(0); Real wy(0);
                        Real uz(0); Real vz(0); Real wz(0);

                        if (flag_fab(i,j,k).isSingleValued())
                        {
                            // Need to check if there are covered cells in neighbours --
                            // -- if so, use one-sided difference computation (but still quadratic)
                            if (flag_fab(i+1,j,k).isCovered())
                            {
                                // Covered cell to the right, go fish left
                                ux = - (c0 * ccvel_fab(i,j,k,0)
                                    + c1 * ccvel_fab(i-1,j,k,0)
                                    + c2 * ccvel_fab(i-2,j,k,0)) * idx;
                                vx = - (c0 * ccvel_fab(i,j,k,1)
                                    + c1 * ccvel_fab(i-1,j,k,1)
                                    + c2 * ccvel_fab(i-2,j,k,1)) * idx;
                                wx = - (c0 * ccvel_fab(i,j,k,2)
                                    + c1 * ccvel_fab(i-1,j,k,2)
                                    + c2 * ccvel_fab(i-2,j,k,2)) * idx;
                            }
                            else if (flag_fab(i-1,j,k).isCovered())
                            {
                                // Covered cell to the left, go fish right
                                ux = (c0 * ccvel_fab(i,j,k,0)
                                    + c1 * ccvel_fab(i+1,j,k,0)
                                    + c2 * ccvel_fab(i+2,j,k,0)) * idx;
                                vx = (c0 * ccvel_fab(i,j,k,1)
                                    + c1 * ccvel_fab(i+1,j,k,1)
                                    + c2 * ccvel_fab(i+2,j,k,1)) * idx;
                                wx = (c0 * ccvel_fab(i,j,k,2)
                                    + c1 * ccvel_fab(i+1,j,k,2)
                                    + c2 * ccvel_fab(i+2,j,k,2)) * idx;
                            }
                            else
                            {
                               // No covered cells right or left, use standard stencil
                                ux = 0.5 * (ccvel_fab(i+1,j,k,0) - ccvel_fab(i-1,j,k,0)) * idx;
                                vx = 0.5 * (ccvel_fab(i+1,j,k,1) - ccvel_fab(i-1,j,k,1)) * idx;
                                wx = 0.5 * (ccvel_fab(i+1,j,k,2) - ccvel_fab(i-1,j,k,2)) * idx;
                            }
                            // Do the same in y-direction
                            if (flag_fab(i,j+1,k).isCovered())
                            {
                                uy = - (c0 * ccvel_fab(i,j,k,0)
                                    + c1 * ccvel_fab(i,j-1,k,0)
                                    + c2 * ccvel_fab(i,j-2,k,0)) * idy;
                                vy = - (c0 * ccvel_fab(i,j,k,1)
                                    + c1 * ccvel_fab(i,j-1,k,1)
                                    + c2 * ccvel_fab(i,j-2,k,1)) * idy;
                                wy = - (c0 * ccvel_fab(i,j,k,2)
                                    + c1 * ccvel_fab(i,j-1,k,2)
                                    + c2 * ccvel_fab(i,j-2,k,2)) * idy;
                            }
                            else if (flag_fab(i,j-1,k).isCovered())
                            {
                                uy = (c0 * ccvel_fab(i,j,k,0)
                                    + c1 * ccvel_fab(i,j+1,k,0)
                                    + c2 * ccvel_fab(i,j+2,k,0)) * idy;
                                vy = (c0 * ccvel_fab(i,j,k,1)
                                    + c1 * ccvel_fab(i,j+1,k,1)
                                    + c2 * ccvel_fab(i,j+2,k,1)) * idy;
                                wy = (c0 * ccvel_fab(i,j,k,2)
                                    + c1 * ccvel_fab(i,j+1,k,2)
                                    + c2 * ccvel_fab(i,j+2,k,2)) * idy;
                            }
                            else
                            {
                                uy = 0.5 * (ccvel_fab(i,j+1,k,0) - ccvel_fab(i,j-1,k,0)) * idy;
                                vy = 0.5 * (ccvel_fab(i,j+1,k,1) - ccvel_fab(i,j-1,k,1)) * idy;
                                wy = 0.5 * (ccvel_fab(i,j+1,k,2) - ccvel_fab(i,j-1,k,2)) * idy;
                            }

                            // Do the same in z-direction
                            if (flag_fab(i,j,k+1).isCovered())
                            {
                                uz = - (c0 * ccvel_fab(i,j,k,0)
                                    + c1 * ccvel_fab(i,j,k-1,0)
                                    + c2 * ccvel_fab(i,j,k-2,0)) * idz;
                                vz = - (c0 * ccvel_fab(i,j,k,1)
                                    + c1 * ccvel_fab(i,j,k-1,1)
                                    + c2 * ccvel_fab(i,j,k-2,1)) * idz;
                                wz = - (c0 * ccvel_fab(i,j,k,2)
                                    + c1 * ccvel_fab(i,j,k-1,2)
                                    + c2 * ccvel_fab(i,j,k-2,2)) * idz;
                            }
                            else if (flag_fab(i,j,k-1).isCovered())
                            {
                                uz = (c0 * ccvel_fab(i,j,k,0)
                                    + c1 * ccvel_fab(i,j,k+1,0)
                                    + c2 * ccvel_fab(i,j,k+2,0)) * idz;
                                vz = (c0 * ccvel_fab(i,j,k,1)
                                    + c1 * ccvel_fab(i,j,k+1,1)
                                    + c2 * ccvel_fab(i,j,k+2,1)) * idz;
                                wz = (c0 * ccvel_fab(i,j,k,2)
                                    + c1 * ccvel_fab(i,j,k+1,2)
                                    + c2 * ccvel_fab(i,j,k+2,2)) * idz;
                            }
                            else
                            {
                                uz = 0.5 * (ccvel_fab(i,j,k+1,0) - ccvel_fab(i,j,k-1,0)) * idz;
                                vz = 0.5 * (ccvel_fab(i,j,k+1,1) - ccvel_fab(i,j,k-1,1)) * idz;
                                wz = 0.5 * (ccvel_fab(i,j,k+1,2) - ccvel_fab(i,j,k-1,2)) * idz;
                            }
                        }
                        else
                        {
                            ux = 0.5 * (ccvel_fab(i+1,j,k,0) - ccvel_fab(i-1,j,k,0)) * idx;
                            vx = 0.5 * (ccvel_fab(i+1,j,k,1) - ccvel_fab(i-1,j,k,1)) * idx;
                            wx = 0.5 * (ccvel_fab(i+1,j,k,2) - ccvel_fab(i-1,j,k,2)) * idx;

                            uy = 0.5 * (ccvel_fab(i,j+1,k,0) - ccvel_fab(i,j-1,k,0)) * idy;
                            vy = 0.5 * (ccvel_fab(i,j+1,k,1) - ccvel_fab(i,j-1,k,1)) * idy;
                            wy = 0.5 * (ccvel_fab(i,j+1,k,2) - ccvel_fab(i,j-1,k,2)) * idy;

                            uz = 0.5 * (ccvel_fab(i,j,k+1,0) - ccvel_fab(i,j,k-1,0)) * idz;
                            vz = 0.5 * (ccvel_fab(i,j,k+1,1) - ccvel_fab(i,j,k-1,1)) * idz;
                            wz = 0.5 * (ccvel_fab(i,j,k+1,2) - ccvel_fab(i,j,k-1,2)) * idz;
                        }
                        sr_fab(i,j,k) = sqrt(2.0 * pow(ux, 2) + 2.0 * pow(vy, 2) + 2.0 * pow(wz, 2)
                                + pow(uy + vx, 2) + pow(vz + wy, 2) + pow(wx + uz, 2));
                    }
                });

                Gpu::synchronize();
            }
        }
    }
}


Real incflo::ComputeKineticEnergy()
{
    BL_PROFILE("incflo::ComputeKineticEnergy");

    // integrated total Kinetic energy
    Real KE = 0.0;

    for(int lev = 0; lev <= finest_level; lev++)
    {
        Real cell_vol = geom[lev].CellSize()[0]*geom[lev].CellSize()[1]*geom[lev].CellSize()[2];

        KE += amrex::ReduceSum(*density[lev],*vel[lev],*level_mask[lev],0,
        [=] AMREX_GPU_HOST_DEVICE (Box const& bx,
                                   Array4<Real const> const& den_arr,
                                   Array4<Real const> const& vel_arr,
                                   Array4<int const>  const& mask_arr) -> Real
        {
            Real KE_Fab = 0.0;

            amrex::Loop(bx, [=,&KE_Fab] (int i, int j, int k) noexcept
            {
                KE_Fab += cell_vol*mask_arr(i,j,k)*den_arr(i,j,k)*( vel_arr(i,j,k,0)*vel_arr(i,j,k,0)
                                                                   +vel_arr(i,j,k,1)*vel_arr(i,j,k,1)
                                                                   +vel_arr(i,j,k,2)*vel_arr(i,j,k,2));

            });
            return KE_Fab;
            
        });
    }
    
    // total volume of grid on level 0
    Real total_vol = geom[0].ProbDomain().volume();
    
    KE *= 0.5/total_vol/ro_0;

    ParallelDescriptor::ReduceRealSum(KE);

    return KE;
    
}

void incflo::ComputeVorticity(Real time_in)
{
	BL_PROFILE("incflo::ComputeVorticity");

    for(int lev = 0; lev <= finest_level; lev++)
    {
        Box domain(geom[lev].Domain());
        Real idx = 1.0 / geom[lev].CellSize()[0];
        Real idy = 1.0 / geom[lev].CellSize()[1];
        Real idz = 1.0 / geom[lev].CellSize()[2];

        // State with ghost cells
        MultiFab Sborder(grids[lev], dmap[lev], vel[lev]->nComp(), nghost,
                         MFInfo(), *ebfactory[lev]);
        FillPatchVel(lev, time_in, Sborder);

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for(MFIter mfi(Sborder, TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            // Tilebox
            Box bx = mfi.tilebox();

            // This is to check efficiently if this tile contains any eb stuff
            const EBFArrayBox& vel_fab = static_cast<EBFArrayBox const&>(Sborder[mfi]);
            const EBCellFlagFab& flags = vel_fab.getEBCellFlagFab();

            // Cell-centered velocity
            const auto& ccvel_fab = Sborder.array(mfi);

            // Cell-centered strain-rate magnitude
            const auto& vort_fab = vort[lev]->array(mfi);

            if (flags.getType(amrex::grow(bx, 0)) == FabType::covered)
            {
                (*vort[lev])[mfi].setVal(1.2345e200, bx);
            }
            else if(flags.getType(amrex::grow(bx, 1)) == FabType::regular)
            {
                // No cut cells in tile + 1-cell witdh halo -> use non-eb routine
                AMREX_FOR_3D(bx, i, j, k,
                {
                    Real vx = 0.5 * (ccvel_fab(i+1,j,k,1) - ccvel_fab(i-1,j,k,1)) * idx;
                    Real wx = 0.5 * (ccvel_fab(i+1,j,k,2) - ccvel_fab(i-1,j,k,2)) * idx;

                    Real uy = 0.5 * (ccvel_fab(i,j+1,k,0) - ccvel_fab(i,j-1,k,0)) * idy;
                    Real wy = 0.5 * (ccvel_fab(i,j+1,k,2) - ccvel_fab(i,j-1,k,2)) * idy;

                    Real uz = 0.5 * (ccvel_fab(i,j,k+1,0) - ccvel_fab(i,j,k-1,0)) * idz;
                    Real vz = 0.5 * (ccvel_fab(i,j,k+1,1) - ccvel_fab(i,j,k-1,1)) * idz;

                    vort_fab(i,j,k) = sqrt( pow(wy - vz, 2) + pow(uz - wx, 2) + pow(vx - uy, 2));
                });
            }
            else
            {
                // Cut cells present -> use EB routine!
                const auto& flag_fab = flags.array();
                Real c0 = -1.5;
                Real c1 = 2.0;
                Real c2 = -0.5;

                AMREX_FOR_3D(bx, i, j, k,
                {
                    if (flag_fab(i,j,k).isCovered())
                    {
                        // Don't compute strainrate in cut cells
                        vort_fab(i,j,k) = 1.2345e200;
                    }
                    else
                    {
                        Real vx(0); Real wx(0);
                        Real uy(0); Real wy(0);
                        Real uz(0); Real vz(0);

                        if (flag_fab(i,j,k).isSingleValued())
                        {
                            // Need to check if there are covered cells in neighbours --
                            // -- if so, use one-sided difference computation (but still quadratic)
                            if (flag_fab(i+1,j,k).isCovered())
                            {
                                // Covered cell to the right, go fish left
                                vx = - (c0 * ccvel_fab(i,j,k,1)
                                    + c1 * ccvel_fab(i-1,j,k,1)
                                    + c2 * ccvel_fab(i-2,j,k,1)) * idx;
                                wx = - (c0 * ccvel_fab(i,j,k,2)
                                    + c1 * ccvel_fab(i-1,j,k,2)
                                    + c2 * ccvel_fab(i-2,j,k,2)) * idx;
                            }
                            else if (flag_fab(i-1,j,k).isCovered())
                            {
                                // Covered cell to the left, go fish right
                                vx = (c0 * ccvel_fab(i,j,k,1)
                                    + c1 * ccvel_fab(i+1,j,k,1)
                                    + c2 * ccvel_fab(i+2,j,k,1)) * idx;
                                wx = (c0 * ccvel_fab(i,j,k,2)
                                    + c1 * ccvel_fab(i+1,j,k,2)
                                    + c2 * ccvel_fab(i+2,j,k,2)) * idx;
                            }
                            else
                            {
                                // No covered cells right or left, use standard stencil
                                vx = 0.5 * (ccvel_fab(i+1,j,k,1) - ccvel_fab(i-1,j,k,1)) * idx;
                                wx = 0.5 * (ccvel_fab(i+1,j,k,2) - ccvel_fab(i-1,j,k,2)) * idx;
                            }
                            // Do the same in y-direction
                            if (flag_fab(i,j+1,k).isCovered())
                            {
                                uy = - (c0 * ccvel_fab(i,j,k,0)
                                    + c1 * ccvel_fab(i,j-1,k,0)
                                    + c2 * ccvel_fab(i,j-2,k,0)) * idy;
                                wy = - (c0 * ccvel_fab(i,j,k,2)
                                    + c1 * ccvel_fab(i,j-1,k,2)
                                    + c2 * ccvel_fab(i,j-2,k,2)) * idy;
                            }
                            else if (flag_fab(i,j-1,k).isCovered())
                            {
                                uy = (c0 * ccvel_fab(i,j,k,0)
                                    + c1 * ccvel_fab(i,j+1,k,0)
                                    + c2 * ccvel_fab(i,j+2,k,0)) * idy;
                                wy = (c0 * ccvel_fab(i,j,k,2)
                                    + c1 * ccvel_fab(i,j+1,k,2)
                                    + c2 * ccvel_fab(i,j+2,k,2)) * idy;
                            }
                            else
                            {
                                uy = 0.5 * (ccvel_fab(i,j+1,k,0) - ccvel_fab(i,j-1,k,0)) * idy;
                                wy = 0.5 * (ccvel_fab(i,j+1,k,2) - ccvel_fab(i,j-1,k,2)) * idy;
                            }

                            // Do the same in z-direction
                            if (flag_fab(i,j,k+1).isCovered())
                            {
                                uz = - (c0 * ccvel_fab(i,j,k,0)
                                    + c1 * ccvel_fab(i,j,k-1,0)
                                    + c2 * ccvel_fab(i,j,k-2,0)) * idz;
                                vz = - (c0 * ccvel_fab(i,j,k,1)
                                    + c1 * ccvel_fab(i,j,k-1,1)
                                    + c2 * ccvel_fab(i,j,k-2,1)) * idz;
                            }
                            else if (flag_fab(i,j,k-1).isCovered())
                            {
                                uz = (c0 * ccvel_fab(i,j,k,0)
                                    + c1 * ccvel_fab(i,j,k+1,0)
                                    + c2 * ccvel_fab(i,j,k+2,0)) * idz;
                                vz = (c0 * ccvel_fab(i,j,k,1)
                                    + c1 * ccvel_fab(i,j,k+1,1)
                                    + c2 * ccvel_fab(i,j,k+2,1)) * idz;
                            }
                            else
                            {
                                uz = 0.5 * (ccvel_fab(i,j,k+1,0) - ccvel_fab(i,j,k-1,0)) * idz;
                                vz = 0.5 * (ccvel_fab(i,j,k+1,1) - ccvel_fab(i,j,k-1,1)) * idz;
                            }
                        }
                        else
                        {
                            vx = 0.5 * (ccvel_fab(i+1,j,k,1) - ccvel_fab(i-1,j,k,1)) * idx;
                            wx = 0.5 * (ccvel_fab(i+1,j,k,2) - ccvel_fab(i-1,j,k,2)) * idx;

                            uy = 0.5 * (ccvel_fab(i,j+1,k,0) - ccvel_fab(i,j-1,k,0)) * idy;
                            wy = 0.5 * (ccvel_fab(i,j+1,k,2) - ccvel_fab(i,j-1,k,2)) * idy;

                            uz = 0.5 * (ccvel_fab(i,j,k+1,0) - ccvel_fab(i,j,k-1,0)) * idz;
                            vz = 0.5 * (ccvel_fab(i,j,k+1,1) - ccvel_fab(i,j,k-1,1)) * idz;
                        }
                        vort_fab(i,j,k) = sqrt(pow(wy-vz,2) + pow(uz-wx,2) + pow(vx-uy,2));
                    }
                });

                Gpu::synchronize();
            } // Cut cells
        } // MFIter
    } // Loop over levels
}

void incflo::ComputeDrag()
{
	BL_PROFILE("incflo::ComputeDrag");

    // Coefficients for one-sided difference estimation
    Real c0 = -1.5;
    Real c1 = 2.0;
    Real c2 = -0.5;

    for(int lev = 0; lev <= finest_level; lev++)
    {
        Box domain(geom[lev].Domain());
        Real dx = geom[lev].CellSize()[0];

        // Get EB geometric info
        const amrex::MultiCutFab* bndryarea;
        const amrex::MultiCutFab* bndrynorm;
        bndryarea = &(ebfactory[lev]->getBndryArea());
        bndrynorm = &(ebfactory[lev]->getBndryNormal());

#ifdef _OPENMP
#pragma omp parallel for reduction(+:drag) if (Gpu::notInLaunchRegion())
#endif
        for(MFIter mfi(*vel[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            // Tilebox
            Box bx = mfi.tilebox();

            // This is to check efficiently if this tile contains any eb stuff
            const EBFArrayBox& vel_fab = static_cast<EBFArrayBox const&>((*vel[lev])[mfi]);
            const EBCellFlagFab& flags = vel_fab.getEBCellFlagFab();

            if (flags.getType(bx) == FabType::singlevalued)
            {
                const auto& drag_arr = drag[lev]->array(mfi);
                const auto& vel_arr = vel[lev]->array(mfi);
                const auto& eta_arr = eta[lev]->array(mfi);
                const auto& p_arr = p[lev]->array(mfi);
                const auto& bndryarea_arr = bndryarea->array(mfi);
                const auto& bndrynorm_arr = bndrynorm->array(mfi);
                const auto& flag_arr = flags.array();

                for(int i = bx.smallEnd(0); i <= bx.bigEnd(0); i++)
                for(int j = bx.smallEnd(1); j <= bx.bigEnd(1); j++)
                for(int k = bx.smallEnd(2); k <= bx.bigEnd(2); k++)
                {
                    if(flag_arr(i,j,k).isSingleValued())
                    {
                        Real area = bndryarea_arr(i,j,k);
                        Real nx = bndrynorm_arr(i,j,k,0);
                        Real ny = bndrynorm_arr(i,j,k,1);
                        Real nz = bndrynorm_arr(i,j,k,2);

                        Real uz, vz, wx, wy, wz;

                        if(flag_arr(i,j,k+1).isCovered())
                        {
                            uz = - (c0 * vel_arr(i,j,k,0) + c1 * vel_arr(i,j,k-1,0) + c2 * vel_arr(i,j,k-2,0)) / dx;
                            vz = - (c0 * vel_arr(i,j,k,1) + c1 * vel_arr(i,j,k-1,1) + c2 * vel_arr(i,j,k-2,1)) / dx;
                            wz = - (c0 * vel_arr(i,j,k,2) + c1 * vel_arr(i,j,k-1,2) + c2 * vel_arr(i,j,k-2,2)) / dx;
                        }
                        else if(flag_arr(i,j,k-1).isCovered())
                        {
                            uz = (c0 * vel_arr(i,j,k,0) + c1 * vel_arr(i,j,k+1,0) + c2 * vel_arr(i,j,k+2,0)) / dx;
                            vz = (c0 * vel_arr(i,j,k,1) + c1 * vel_arr(i,j,k+1,1) + c2 * vel_arr(i,j,k+2,1)) / dx;
                            wz = (c0 * vel_arr(i,j,k,2) + c1 * vel_arr(i,j,k+1,2) + c2 * vel_arr(i,j,k+2,2)) / dx;
                        }
                        else
                        {
                            uz = 0.5 * (vel_arr(i,j,k+1,0) - vel_arr(i,j,k-1,0)) / dx;
                            vz = 0.5 * (vel_arr(i,j,k+1,1) - vel_arr(i,j,k-1,1)) / dx;
                            wz = 0.5 * (vel_arr(i,j,k+1,2) - vel_arr(i,j,k-1,2)) / dx;
                        }

                        if(flag_arr(i,j+1,k).isCovered())
                        {
                            wy = - (c0 * vel_arr(i,j,k,2) + c1 * vel_arr(i,j-1,k,2) + c2 * vel_arr(i,j-2,k,2)) / dx;
                        }
                        else if(flag_arr(i,j-1,k).isCovered())
                        {
                            wy = (c0 * vel_arr(i,j,k,2) + c1 * vel_arr(i,j+1,k,2) + c2 * vel_arr(i,j+2,k,2)) / dx;
                        }
                        else
                        {
                            wy = 0.5 * (vel_arr(i,j+1,k,2) - vel_arr(i,j-1,k,2)) / dx;
                        }

                        if(flag_arr(i+1,j,k).isCovered())
                        {
                            wx = - (c0 * vel_arr(i,j,k,2) + c1 * vel_arr(i-1,j,k,2) + c2 * vel_arr(i-2,j,k,2)) / dx;
                        }
                        else if(flag_arr(i-1,j,k).isCovered())
                        {
                            wx = (c0 * vel_arr(i,j,k,2) + c1 * vel_arr(i+1,j,k,2) + c2 * vel_arr(i+2,j,k,2)) / dx;
                        }
                        else
                        {
                            wx = 0.5 * (vel_arr(i+1,j,k,2) - vel_arr(i-1,j,k,2)) / dx;
                        }

                        Real p_contrib = p_arr(i,j,k) * nz;
                        Real tau_contrib = - eta_arr(i,j,k) * ( (uz + wx) * nx + (vz + wy) * ny + (wz + wz) * nz );

                        // TODO: Get values on EB centroid,
                        //       not the default CC and nodal values
                        drag_arr(i,j,k) = (p_contrib + tau_contrib) * area * dx * dx;
                    }
                    else
                    {
                        drag_arr(i,j,k) = 0.0;
                    }
                }
            }
            else
            {
                (*drag[lev])[mfi].setVal(0.0, bx);
            }
        }
    }
}
