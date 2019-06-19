#include <AMReX_Box.H>
#include <AMReX_EBMultiFabUtil.H>

#include <incflo.H>
#include <derive_F.H>
#include <projection_F.H>

void incflo::UpdateDerivedQuantities()
{
    BL_PROFILE("incflo::UpdateDerivedQuantities()");

    ComputeDivU(cur_time);
    ComputeStrainrate();
    ComputeViscosity();
    ComputeVorticity();
}

void incflo::ComputeDivU(Real time)
{
    int extrap_dir_bcs = 0;
    FillVelocityBC(time, extrap_dir_bcs);

    // Define the operator in order to compute the multi-level divergence
    //
    //        (del dot b sigma grad)) phi
    //
    LPInfo info;
    MLNodeLaplacian matrix(geom, grids, dmap, info, amrex::GetVecOfConstPtrs(ebfactory));

    // Set domain BCs for Poisson's solver
    // The domain BCs refer to level 0 only
    int bc_lo[3], bc_hi[3];
    Box domain(geom[0].Domain());

    set_ppe_bc(bc_lo, bc_hi,
               domain.loVect(), domain.hiVect(),
               &nghost,
               bc_ilo[0]->dataPtr(), bc_ihi[0]->dataPtr(),
               bc_jlo[0]->dataPtr(), bc_jhi[0]->dataPtr(),
               bc_klo[0]->dataPtr(), bc_khi[0]->dataPtr());

    matrix.setDomainBC({(LinOpBCType)bc_lo[0], (LinOpBCType)bc_lo[1], (LinOpBCType)bc_lo[2]},
                       {(LinOpBCType)bc_hi[0], (LinOpBCType)bc_hi[1], (LinOpBCType)bc_hi[2]});

    matrix.compDivergence(GetVecOfPtrs(divu), GetVecOfPtrs(vel)); 
}

void incflo::ComputeStrainrate()
{
    BL_PROFILE("incflo::ComputeStrainrate");

    for(int lev = 0; lev <= finest_level; lev++)
    {
        Box domain(geom[lev].Domain());

        // State with ghost cells
        MultiFab Sborder(grids[lev], dmap[lev], vel[lev]->nComp(), nghost, 
                         MFInfo(), *ebfactory[lev]);
        FillPatchVel(lev, cur_time, Sborder, 0, Sborder.nComp());
    
        // Copy each FAB back from Sborder into the vel array, complete with filled ghost cells
        MultiFab::Copy(*vel[lev], Sborder, 0, 0, vel[lev]->nComp(), vel[lev]->nGrow());

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

            if (flags.getType(bx) == FabType::covered)
            {
                (*strainrate[lev])[mfi].setVal(1.2345e200, bx);
            }
            else
            {
                if(flags.getType(amrex::grow(bx, 0)) == FabType::regular)
                {
                    compute_strainrate(BL_TO_FORTRAN_BOX(bx),
                                       BL_TO_FORTRAN_ANYD((*strainrate[lev])[mfi]),
                                       BL_TO_FORTRAN_ANYD((*vel[lev])[mfi]),
                                       geom[lev].CellSize());
                }
                else
                {
                    compute_strainrate_eb(BL_TO_FORTRAN_BOX(bx),
                                          BL_TO_FORTRAN_ANYD((*strainrate[lev])[mfi]),
                                          BL_TO_FORTRAN_ANYD((*vel[lev])[mfi]),
                                          BL_TO_FORTRAN_ANYD(flags),
                                          geom[lev].CellSize());
                }
            }
        }
    }
}

void incflo::ComputeVorticity()
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
        FillPatchVel(lev, cur_time, Sborder, 0, Sborder.nComp());
    
        // Copy each FAB back from Sborder into the vel array, complete with filled ghost cells
        MultiFab::Copy (*vel[lev], Sborder, 0, 0, vel[lev]->nComp(), vel[lev]->nGrow());

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

            if (flags.getType(bx) == FabType::covered)
            {
                (*vort[lev])[mfi].setVal(1.2345e200, bx);
            }
            else
            {
                if(flags.getType(amrex::grow(bx, 0)) == FabType::regular)
                {
                    const auto& vel_arr = Sborder.array(mfi);
                    const auto& vort_arr = vort[lev]->array(mfi);

                    for(int i = bx.smallEnd(0); i <= bx.bigEnd(0); i++)
                    for(int j = bx.smallEnd(1); j <= bx.bigEnd(1); j++)
                    for(int k = bx.smallEnd(2); k <= bx.bigEnd(2); k++)
                    {
                        Real vx = (vel_arr(i+1, j  , k  , 1) - vel_arr(i-1, j  , k  , 1)) * idx;
                        Real wx = (vel_arr(i+1, j  , k  , 2) - vel_arr(i-1, j  , k  , 2)) * idx;
                        Real uy = (vel_arr(i  , j+1, k  , 0) - vel_arr(i  , j-1, k  , 0)) * idy;
                        Real wy = (vel_arr(i  , j+1, k  , 2) - vel_arr(i  , j-1, k  , 2)) * idy;
                        Real uz = (vel_arr(i  , j  , k+1, 0) - vel_arr(i  , j  , k-1, 0)) * idz;
                        Real vz = (vel_arr(i  , j  , k+1, 1) - vel_arr(i  , j  , k-1, 1)) * idz;
                        
                        // The factor half is included here instead of in each of the above
                        vort_arr(i,j,k) = 0.5 * sqrt(pow(wy - vz, 2) + pow(uz - wx, 2) + pow(vx - uy, 2));
                    }
                }
                else
                {
                    compute_vort_eb(BL_TO_FORTRAN_BOX(bx),
                                    BL_TO_FORTRAN_ANYD((*vort[lev])[mfi]),
                                    BL_TO_FORTRAN_ANYD((*vel[lev])[mfi]),
                                    BL_TO_FORTRAN_ANYD(flags),
                                    geom[lev].CellSize());
                }
            }
        }
    }
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

                        Real ux, vx, wx, uy, uz;

                        if(flag_arr(i+1,j,k).isCovered())
                        {
                            ux = - (c0 * vel_arr(i,j,k,0) + c1 * vel_arr(i-1,j,k,0) + c2 * vel_arr(i-2,j,k,0)) / dx;
                            vx = - (c0 * vel_arr(i,j,k,1) + c1 * vel_arr(i-1,j,k,1) + c2 * vel_arr(i-2,j,k,1)) / dx;
                            wx = - (c0 * vel_arr(i,j,k,2) + c1 * vel_arr(i-1,j,k,2) + c2 * vel_arr(i-2,j,k,2)) / dx;
                        }
                        else if(flag_arr(i-1,j,k).isCovered())
                        {
                            ux = (c0 * vel_arr(i,j,k,0) + c1 * vel_arr(i+1,j,k,0) + c2 * vel_arr(i+2,j,k,0)) / dx;
                            vx = (c0 * vel_arr(i,j,k,1) + c1 * vel_arr(i+1,j,k,1) + c2 * vel_arr(i+2,j,k,1)) / dx;
                            wx = (c0 * vel_arr(i,j,k,2) + c1 * vel_arr(i+1,j,k,2) + c2 * vel_arr(i+2,j,k,2)) / dx;
                        }
                        else
                        {
                            ux = 0.5 * (vel_arr(i+1,j,k,0) - vel_arr(i-1,j,k,0)) / dx;
                            vx = 0.5 * (vel_arr(i+1,j,k,1) - vel_arr(i-1,j,k,1)) / dx;
                            wx = 0.5 * (vel_arr(i+1,j,k,2) - vel_arr(i-1,j,k,2)) / dx;
                        }

                        if(flag_arr(i,j+1,k).isCovered())
                        {
                            uy = - (c0 * vel_arr(i,j,k,0) + c1 * vel_arr(i,j-1,k,0) + c2 * vel_arr(i,j-2,k,0)) / dx;
                        }
                        else if(flag_arr(i,j-1,k).isCovered())
                        {
                            uy = (c0 * vel_arr(i,j,k,0) + c1 * vel_arr(i,j+1,k,0) + c2 * vel_arr(i,j+2,k,0)) / dx;
                        }
                        else
                        {
                            uy = 0.5 * (vel_arr(i,j+1,k,0) - vel_arr(i,j-1,k,0)) / dx;
                        }

                        if(flag_arr(i,j,k+1).isCovered())
                        {
                            uz = - (c0 * vel_arr(i,j,k,0) + c1 * vel_arr(i,j,k-1,0) + c2 * vel_arr(i,j,k-2,0)) / dx;
                        }
                        else if(flag_arr(i,j,k-1).isCovered())
                        {
                            uz = (c0 * vel_arr(i,j,k,0) + c1 * vel_arr(i,j,k+1,0) + c2 * vel_arr(i,j,k+2,0)) / dx;
                        }
                        else
                        {
                            uz = 0.5 * (vel_arr(i,j,k+1,0) - vel_arr(i,j,k-1,0)) / dx;
                        }

                        Real p_contrib = p_arr(i,j,k) * nx;
                        Real tau_contrib = - eta_arr(i,j,k) * ( (ux + ux) * nx + (vx + uy) * ny + (wx + uz) * nz );

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
