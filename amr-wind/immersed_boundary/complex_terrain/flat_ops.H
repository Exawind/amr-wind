#ifndef FLAT_OPS_H
#define FLAT_OPS_H

#include "amr-wind/immersed_boundary/complex_terrain/Flat.H"
#include "amr-wind/immersed_boundary/IBOps.H"
#include "amr-wind/immersed_boundary/IB.H"
#include "amr-wind/immersed_boundary/complex_terrain/complex_terrain_ops.H"

namespace amr_wind {
namespace ib {

namespace ops {

template <>
struct ReadInputsOp<Flat>
{
    void
    operator()(Flat::DataType& data, const ::amr_wind::utils::MultiParser& pp)
    {
        auto& wdata = data.meta();
        auto& info = data.info();

        complex_terrain::read_inputs(wdata, info, pp);
        pp.get("height", wdata.height);
    }
};

template <>
struct InitDataOp<Flat>
{
    void operator()(Flat::DataType& data)
    {
        const auto& wdata = data.meta();

        auto& sim = data.sim();
        // cppcheck-suppress constVariable
        auto& mask_node = sim.repo().get_int_field("mask_node");
        // cppcheck-suppress constVariable
        auto& levelset = sim.repo().get_field("ib_levelset");

        auto nlevels = sim.repo().num_active_levels();
        auto geom = sim.mesh().Geom();

        amrex::Real h = wdata.height;

        for (int lev = 0; lev < nlevels; ++lev) {
            const auto& problo = geom[lev].ProbLoArray();
            const auto& dx = geom[lev].CellSizeArray();

            for (amrex::MFIter mfi(levelset(lev)); mfi.isValid(); ++mfi) {
                const auto& bx = mfi.growntilebox();
                auto phi = levelset(lev).array(mfi);
                amrex::ParallelFor(
                    bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                        amrex::Real z = problo[2] + (k + 0.5) * dx[2];
                        const amrex::Real phi_glob = phi(i, j, k);
                        amrex::Real phi_loc = z - h;
                        phi(i, j, k) = std::min(phi_loc, phi_glob);
                    });

                const auto& nbx = mfi.nodaltilebox();
                auto epsilon_node = mask_node(lev).array(mfi);
                amrex::ParallelFor(
                    nbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                        const amrex::Real z = problo[2] + k * dx[2];
                        if (z <= h) {
                            epsilon_node(i, j, k) = 0;
                        }
                    });
            }
        }
    }
};

} // namespace ops
} // namespace ib
} // namespace amr_wind

#endif /* FLAT_OPS_H */
