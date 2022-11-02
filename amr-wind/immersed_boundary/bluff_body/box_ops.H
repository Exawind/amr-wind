#ifndef BOX_OPS_H
#define BOX_OPS_H

#include "amr-wind/immersed_boundary/bluff_body/Box.H"
#include "amr-wind/immersed_boundary/IBOps.H"
#include "amr-wind/immersed_boundary/IB.H"
#include "amr-wind/immersed_boundary/bluff_body/bluff_body_ops.H"

namespace amr_wind::ib::ops {

template <>
struct ReadInputsOp<Box>
{
    void
    operator()(Box::DataType& data, const ::amr_wind::utils::MultiParser& pp)
    {
        auto& wdata = data.meta();
        auto& info = data.info();

        bluff_body::read_inputs(wdata, info, pp);

        pp.get("center", wdata.center_loc);
        pp.get("length", wdata.length);
        pp.get("width", wdata.width);
        pp.get("height", wdata.height);

        amrex::Real search_radius =
            2.0 * std::max(wdata.length, std::max(wdata.width, wdata.height));

        // clang-format off
        const auto& origin=wdata.center_loc;
        info.bound_box = amrex::RealBox(
            origin[0] - search_radius,
            origin[1] - search_radius,
            origin[2] - search_radius,
            origin[0] + search_radius,
            origin[1] + search_radius,
            origin[2] + search_radius);
        // clang-format on
    }
};

template <>
struct InitDataOp<Box>
{
    void operator()(Box::DataType& data)
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
                auto epsilon_node = mask_node(lev).array(mfi);
                auto phi = levelset(lev).array(mfi);

                const amrex::Real x0 = wdata.center_loc[0];
                const amrex::Real y0 = wdata.center_loc[1];
                const amrex::Real z0 = wdata.center_loc[2];
                const amrex::Real l = wdata.length;
                const amrex::Real w = wdata.width;
                const amrex::Real h = wdata.height;
                amrex::ParallelFor(
                    bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                        const amrex::Real x = problo[0] + (i + 0.5) * dx[0];
                        const amrex::Real y = problo[1] + (j + 0.5) * dx[1];
                        const amrex::Real z = problo[2] + (k + 0.5) * dx[2];

                        if ((std::abs(x - x0) <= 0.5 * l) &&
                            (std::abs(y - y0) <= 0.5 * w) &&
                            (std::abs(z - z0) <= 0.5 * h)) {
                            phi(i, j, k) = std::max(
                                std::abs(x - x0) - 0.5 * l,
                                std::max(
                                    std::abs(y - y0) - 0.5 * w,
                                    std::abs(z - z0) - 0.5 * h));
                        } else {
                            phi(i, j, k) = std::min(
                                std::abs(std::abs(x - x0) - 0.5 * l),
                                std::min(
                                    std::abs(std::abs(y - y0) - 0.5 * w),
                                    std::abs(std::abs(z - z0) - 0.5 * h)));
                        }
                    });
                const auto& nbx = mfi.nodaltilebox();
                amrex::ParallelFor(
                    nbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                        const amrex::Real x = problo[0] + i * dx[0];
                        const amrex::Real y = problo[1] + j * dx[1];
                        const amrex::Real z = problo[2] + k * dx[2];

                        if ((std::abs(x - x0) <= 0.5 * l + 0.5 * dx[0]) &&
                            (std::abs(y - y0) <= 0.5 * w + 0.5 * dx[1]) &&
                            (std::abs(z - z0) <= 0.5 * h + 0.5 * dx[2])) {
                            epsilon_node(i, j, k) = 0;
                        }
                    });
            }
        }
    }
};

} // namespace amr_wind::ib::ops

#endif /* BOX_OPS_H */
