#ifndef SPHERE_OPS_H
#define SPHERE_OPS_H

#include "amr-wind/immersed_boundary/bluff_body/Sphere.H"
#include "amr-wind/immersed_boundary/IBOps.H"
#include "amr-wind/immersed_boundary/IB.H"
#include "amr-wind/immersed_boundary/bluff_body/bluff_body_ops.H"

namespace amr_wind {
namespace ib {

namespace ops {

template <>
struct ReadInputsOp<Sphere>
{
    void
    operator()(Sphere::DataType& data, const ::amr_wind::utils::MultiParser& pp)
    {
        auto& wdata = data.meta();
        auto& info = data.info();

        bluff_body::read_inputs(wdata, info, pp);

        pp.get("center", wdata.center_loc);
        pp.get("radius", wdata.radius);

        amrex::Real search_radius = 3.0 * wdata.radius;
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
struct InitDataOp<Sphere>
{
    void operator()(Sphere::DataType& data)
    {
        auto& wdata = data.meta();

        auto& sim = data.sim();
        auto& mask_node = sim.repo().get_int_field("mask_node");
        auto& levelset = sim.repo().get_field("ib_levelset");

        auto nlevels = sim.repo().num_active_levels();
        auto geom = sim.mesh().Geom();

        for (int lev = 0; lev < nlevels; ++lev) {
            const auto& problo = geom[lev].ProbLoArray();
            const auto& dx = geom[lev].CellSizeArray();

            for (amrex::MFIter mfi(levelset(lev)); mfi.isValid(); ++mfi) {
                const auto& bx = mfi.growntilebox();
                auto phi = levelset(lev).array(mfi);

                amrex::Real x0 = wdata.center_loc[0];
                amrex::Real y0 = wdata.center_loc[1];
                amrex::Real z0 = wdata.center_loc[2];
                amrex::Real R = wdata.radius;
                amrex::ParallelFor(
                    bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                        const amrex::Real x = problo[0] + (i + 0.5) * dx[0];
                        const amrex::Real y = problo[1] + (j + 0.5) * dx[1];
                        const amrex::Real z = problo[2] + (k + 0.5) * dx[2];
                        amrex::Real phi_glob = phi(i, j, k);
                        amrex::Real r = std::sqrt(
                            (x - x0) * (x - x0) + (y - y0) * (y - y0) +
                            (z - z0) * (z - z0));

                        amrex::Real phi_loc = r - R;
                        phi(i, j, k) = std::min(phi_loc, phi_glob);
                    });

                const auto& nbx = mfi.nodaltilebox();
                auto epsilon_node = mask_node(lev).array(mfi);
                amrex::ParallelFor(
                    nbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                        const amrex::Real x = problo[0] + i * dx[0];
                        const amrex::Real y = problo[1] + j * dx[1];
                        const amrex::Real z = problo[2] + k * dx[2];

                        amrex::Real r = std::sqrt(
                            (x - x0) * (x - x0) + (y - y0) * (y - y0) +
                            (z - z0) * (z - z0));
                        if (r <= R) {
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
