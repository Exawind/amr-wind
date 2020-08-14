#include <AMReX_Box.H>

#include "amr-wind/incflo.H"
#include <AMReX_NodalProjector.H>

using namespace amrex;

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
