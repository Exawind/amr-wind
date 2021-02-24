#include <AMReX_Box.H>

#include "amr-wind/incflo.H"
#include <AMReX_NodalProjector.H>

using namespace amrex;

Real incflo::ComputeKineticEnergy() const
{
    BL_PROFILE("amr-wind::incflo::ComputeKineticEnergy");

    // integrated total Kinetic energy
    Real KE = 0.0;

    for (int lev = 0; lev <= finest_level; lev++) {

        iMultiFab level_mask;
        if (lev < finest_level) {
            level_mask = makeFineMask(
                grids[lev], dmap[lev], grids[lev + 1], IntVect(2), 1, 0);
        } else {
            level_mask.define(grids[lev], dmap[lev], 1, 0, MFInfo());
            level_mask.setVal(1);
        }

        const Real cell_vol = geom[lev].CellSize()[0] *
                              geom[lev].CellSize()[1] * geom[lev].CellSize()[2];

        KE += amrex::ReduceSum(
            density()(lev), velocity()(lev), level_mask, 0,
            [=] AMREX_GPU_HOST_DEVICE(
                Box const& bx, Array4<Real const> const& den_arr,
                Array4<Real const> const& vel_arr,
                Array4<int const> const& mask_arr) -> Real {
                Real KE_Fab = 0.0;

                amrex::Loop(bx, [=, &KE_Fab](int i, int j, int k) noexcept {
                    KE_Fab += cell_vol * mask_arr(i, j, k) * den_arr(i, j, k) *
                              (vel_arr(i, j, k, 0) * vel_arr(i, j, k, 0) +
                               vel_arr(i, j, k, 1) * vel_arr(i, j, k, 1) +
                               vel_arr(i, j, k, 2) * vel_arr(i, j, k, 2));
                });
                return KE_Fab;
            });
    }

    // total volume of grid on level 0
    const Real total_vol = geom[0].ProbDomain().volume();

    KE *= 0.5 / total_vol;

    ParallelDescriptor::ReduceRealSum(KE);

    return KE;
}
