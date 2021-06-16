#ifndef CYLINDER_OPS_H
#define CYLINDER_OPS_H

#include "amr-wind/immersed_boundary/bluff_body/Cylinder.H"
#include "amr-wind/immersed_boundary/IBOps.H"
#include "amr-wind/immersed_boundary/IB.H"
#include "amr-wind/immersed_boundary/bluff_body/bluff_body_ops.H"

namespace amr_wind {
namespace ib {

namespace ops {

template <>
struct ReadInputsOp<Cylinder>
{
    void operator()(
        Cylinder::DataType& data, const ::amr_wind::utils::MultiParser& pp)
    {

        auto& wdata = data.meta();
        auto& info = data.info();

        bluff_body::read_inputs(wdata, info, pp);

        pp.get("center", wdata.center_loc);
        pp.get("radius", wdata.radius);
        pp.get("length", wdata.length);

        amrex::Real search_radius = 2.0 * wdata.radius;
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
struct InitDataOp<Cylinder>
{
    void operator()(Cylinder::DataType& data)
    {
        auto& wdata = data.meta();

        auto& sim = data.sim();
        auto& mask_cell = sim.repo().get_int_field("mask_cell");
        auto& mask_node = sim.repo().get_int_field("mask_node");

        auto nlevels = sim.repo().num_active_levels();
        auto geom = sim.mesh().Geom();

        for (int lev = 0; lev < nlevels; ++lev) {
            const auto& problo = geom[lev].ProbLoArray();
            const auto& dx = geom[lev].CellSizeArray();
            for (amrex::MFIter mfi(mask_cell(lev)); mfi.isValid(); ++mfi) {
                auto epsilon_cell = mask_cell(lev).array(mfi);
                auto epsilon_node = mask_node(lev).array(mfi);
                amrex::Real x0 = wdata.center_loc[0];
                amrex::Real y0 = wdata.center_loc[1];
                amrex::Real z0 = wdata.center_loc[2];
                amrex::Real l = wdata.length;

                const auto& bx = mfi.growntilebox();
                amrex::ParallelFor(
                    bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                        const amrex::Real x = problo[0] + (i + 0.5) * dx[0];
                        const amrex::Real y = problo[1] + (j + 0.5) * dx[1];
                        const amrex::Real z = problo[2] + (k + 0.5) * dx[2];

                        amrex::Real r = std::sqrt(
                            (x - x0) * (x - x0) + (y - y0) * (y - y0));
                        if (r <= wdata.radius &&
                            (std::abs(z - z0) <= 0.5 * l)) {
                            epsilon_cell(i, j, k) = 0;
                        }
                    });
                const auto& nbx = mfi.nodaltilebox();
                amrex::ParallelFor(
                    nbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                        const amrex::Real x = problo[0] + i * dx[0];
                        const amrex::Real y = problo[1] + j * dx[1];
                        const amrex::Real z = problo[2] + k * dx[2];

                        amrex::Real r = std::sqrt(
                            (x - x0) * (x - x0) + (y - y0) * (y - y0));
                        if (r <= wdata.radius &&
                            (std::abs(z - z0) <= 0.5 * l)) {
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
