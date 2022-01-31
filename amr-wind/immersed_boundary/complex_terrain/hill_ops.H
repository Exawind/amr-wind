#ifndef HILL_OPS_H
#define HILL_OPS_H

#include "amr-wind/immersed_boundary/complex_terrain/Hill.H"
#include "amr-wind/immersed_boundary/IBOps.H"
#include "amr-wind/immersed_boundary/IB.H"
#include "amr-wind/immersed_boundary/complex_terrain/complex_terrain_ops.H"

namespace amr_wind {
namespace ib {

namespace ops {

template <>
struct ReadInputsOp<Hill>
{
    void
    operator()(Hill::DataType& data, const ::amr_wind::utils::MultiParser& pp)
    {
        auto& wdata = data.meta();
        auto& info = data.info();

        complex_terrain::read_inputs(wdata, info, pp);

        pp.get("center", wdata.center_loc);
        pp.get("height", wdata.height);
        pp.get("half_width", wdata.half_width);
    }
};

template <>
struct InitDataOp<Hill>
{
    void operator()(Hill::DataType& data)
    {
        const auto& wdata = data.meta();

        auto& sim = data.sim();
        // cppcheck-suppress constVariable
        auto& mask_node = sim.repo().get_int_field("mask_node");
        // cppcheck-suppress constVariable
        auto& levelset = sim.repo().get_field("ib_levelset");

        auto nlevels = sim.repo().num_active_levels();
        auto geom = sim.mesh().Geom();

        for (int lev = 0; lev < nlevels; ++lev) {
            const auto& problo = geom[lev].ProbLoArray();
            const auto& dx = geom[lev].CellSizeArray();

            for (amrex::MFIter mfi(levelset(lev)); mfi.isValid(); ++mfi) {
                const auto& bx = mfi.growntilebox();
                auto phi = levelset(lev).array(mfi);

                const amrex::Real x0 = wdata.center_loc[0];
                const amrex::Real y0 = wdata.center_loc[1];
                const amrex::Real z0 = wdata.center_loc[2];
                const amrex::Real h = wdata.height;
                const amrex::Real b = wdata.half_width;
                amrex::ParallelFor(
                    bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                        const amrex::Real x = problo[0] + (i + 0.5) * dx[0];
                        const amrex::Real y = problo[1] + (j + 0.5) * dx[1];
                        const amrex::Real z = problo[2] + (k + 0.5) * dx[2];
                        const amrex::Real phi_glob = phi(i, j, k);
                        const amrex::Real dist = std::sqrt(
                            (x - x0) * (x - x0) + (y - y0) * (y - y0) +
                            (z - z0) * (z - z0));
                        const amrex::Real hill =
                            h * std::pow(
                                    std::cos(
                                        M_PI *
                                        std::sqrt(
                                            (x - x0) * (x - x0) +
                                            (y - y0) * (y - y0)) /
                                        (2. * b)),
                                    2.);

                        amrex::Real phi_loc = dist - hill;
                        phi(i, j, k) = std::min(phi_loc, phi_glob);
                    });

                const auto& nbx = mfi.nodaltilebox();
                auto epsilon_node = mask_node(lev).array(mfi);
                amrex::ParallelFor(
                    nbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                        const amrex::Real x = problo[0] + i * dx[0];
                        const amrex::Real y = problo[1] + j * dx[1];
                        const amrex::Real z = problo[2] + k * dx[2];

                        const amrex::Real dist = std::sqrt(
                            (x - x0) * (x - x0) + (y - y0) * (y - y0) +
                            (z - z0) * (z - z0));
                        const amrex::Real hill =
                            h * std::pow(
                                    std::cos(
                                        M_PI *
                                        std::sqrt(
                                            (x - x0) * (x - x0) +
                                            (y - y0) * (y - y0)) /
                                        (2. * b)),
                                    2.);
                        if (dist <= hill) {
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

#endif /* CYLINDER_OPS_H */
