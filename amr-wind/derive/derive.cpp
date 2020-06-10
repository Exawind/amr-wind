#include <AMReX_Box.H>

#include "amr-wind/incflo.H"
#include <AMReX_NodalProjector.H>

using namespace amrex;

void incflo::ComputeDivU(Real /* time_in */)
{
#if 0 //fixme
    BL_PROFILE("amr-wind::incflo::ComputeDivU");

    incflo_set_velocity_bcs(time_in, vel);

    int bc_lo[3], bc_hi[3];
    Box domain(geom[0].Domain());

    set_ppe_bcs(bc_lo, bc_hi,
                domain.loVect(), domain.hiVect(),
                &nghost,
                bc_ilo[0]->dataPtr(), bc_ihi[0]->dataPtr(),
                bc_jlo[0]->dataPtr(), bc_jhi[0]->dataPtr(),
                bc_klo[0]->dataPtr(), bc_khi[0]->dataPtr());

    ppe_lobc = {(LinOpBCType)bc_lo[0], (LinOpBCType)bc_lo[1], (LinOpBCType)bc_lo[2]};
    ppe_hibc = {(LinOpBCType)bc_hi[0], (LinOpBCType)bc_hi[1], (LinOpBCType)bc_hi[2]};

    LPInfo lpinfo;

    //
    // This rebuilds integrals each time linop is created -- must find a better way
    //

    MLNodeLaplacian linop(geom, grids, dmap, lpinfo);

    linop.setDomainBC(ppe_lobc,ppe_hibc);
    linop.compDivergence(GetVecOfPtrs(divu),GetVecOfPtrs(vel));
#endif
}

void incflo::ComputeStrainrate(Real /* time_in */)
{
#if 0
    BL_PROFILE("amr-wind::incflo::ComputeStrainrate");

    for(int lev = 0; lev <= finest_level; lev++)
    {
        Box domain(geom[lev].Domain());
        Real idx = 1.0 / geom[lev].CellSize()[0];
        Real idy = 1.0 / geom[lev].CellSize()[1];
        Real idz = 1.0 / geom[lev].CellSize()[2];

        // State with ghost cells
        MultiFab Sborder(grids[lev], dmap[lev], vel[lev]->nComp(), 1, MFInfo());

        FillPatchVel(lev, time_in, Sborder);

        strainrate[lev]->setVal(1.2345e200);

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(Sborder, TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            // Tilebox
            Box bx = mfi.tilebox();

            // Cell-centered velocity
            const auto& ccvel_fab = Sborder.array(mfi);

            // Cell-centered strain-rate magnitude
            const auto& sr_fab = strainrate[lev]->array(mfi);
            
            amrex::ParallelFor(bx,
              [idx,idy,idz,sr_fab,ccvel_fab] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
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
//fixme use strain rate in derive_K.H
                // Include the factor half here rather than in each of the above
                sr_fab(i,j,k) = sqrt(2.0 * pow(ux, 2) + 2.0 * pow(vy, 2) + 2.0 * pow(wz, 2)
                        + pow(uy + vx, 2) + pow(vz + wy, 2) + pow(wx + uz, 2));
            });
        } // MFIter
    } // lev
#endif
}


Real incflo::ComputeKineticEnergy () const
{
    BL_PROFILE("amr-wind::incflo::ComputeKineticEnergy");

    // integrated total Kinetic energy
    Real KE = 0.0;

    for(int lev = 0; lev <= finest_level; lev++)
    {

        iMultiFab level_mask;
        if(lev < finest_level){
            level_mask = makeFineMask(grids[lev],dmap[lev], grids[lev+1], IntVect(2), 1, 0);
        }else {
            level_mask.define(grids[lev], dmap[lev], 1, 0, MFInfo());
            level_mask.setVal(1);
        }

        const Real cell_vol = geom[lev].CellSize()[0]*geom[lev].CellSize()[1]*geom[lev].CellSize()[2];

        KE += amrex::ReduceSum(density()(lev), velocity()(lev), level_mask, 0,
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
    const Real total_vol = geom[0].ProbDomain().volume();

    KE *= 0.5/total_vol;

    ParallelDescriptor::ReduceRealSum(KE);

    return KE;

}

void incflo::ComputeVorticity (int lev, Real /* t */, MultiFab& vort, MultiFab const& vel)
{
    BL_PROFILE("amr-wind::incflo::ComputeVorticity");

    const Real idx = Geom(lev).InvCellSize(0);
    const Real idy = Geom(lev).InvCellSize(1);
    const Real idz = Geom(lev).InvCellSize(2);

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for(MFIter mfi(vel, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Box bx = mfi.tilebox();
        Array4<Real const> const& ccvel_fab = vel.const_array(mfi);
        Array4<Real> const& vort_fab = vort.array(mfi);
  
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real vx = 0.5 * (ccvel_fab(i+1,j,k,1) - ccvel_fab(i-1,j,k,1)) * idx;
            Real wx = 0.5 * (ccvel_fab(i+1,j,k,2) - ccvel_fab(i-1,j,k,2)) * idx;

            Real uy = 0.5 * (ccvel_fab(i,j+1,k,0) - ccvel_fab(i,j-1,k,0)) * idy;
            Real wy = 0.5 * (ccvel_fab(i,j+1,k,2) - ccvel_fab(i,j-1,k,2)) * idy;

            Real uz = 0.5 * (ccvel_fab(i,j,k+1,0) - ccvel_fab(i,j,k-1,0)) * idz;
            Real vz = 0.5 * (ccvel_fab(i,j,k+1,1) - ccvel_fab(i,j,k-1,1)) * idz;

            vort_fab(i,j,k) = std::sqrt((wy-vz)*(wy-vz) + (uz-wx)*(uz-wx) + (vx-uy)*(vx-uy));
        });
        
    }
}
