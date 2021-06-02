#ifndef SPHERE_OPS_H
#define SPHERE_OPS_H

#include "amr-wind/immersed_boundary/geometry/TaylorGreen.H"
#include "amr-wind/immersed_boundary/IBOps.H"
#include "amr-wind/immersed_boundary/geometry/bluff_body_ops.H"

namespace amr_wind {
namespace ib {

namespace ops {

template <>
struct ReadInputsOp<TaylorGreen>
{
    void operator()(TaylorGreen::DataType& data, const utils::IBParser& pp)
    {
        auto& wdata = data.meta();
        auto& info = data.info();
        pp.get("num_points", wdata.num_pts);
        pp.get("center", wdata.center_loc);
        pp.get("radius", wdata.radius);
        amrex::Real search_radius = 3 * wdata.radius;

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
struct InitDataOp<TaylorGreen>
{
    void operator()(TaylorGreen::DataType& data)
    {
        auto& wdata = data.meta();
        auto& grid = data.grid();
        int npts = wdata.num_pts;
        grid.resize(npts * npts);

        int ip = 0;
        for (int iphi = 0; iphi < npts; ++iphi) {
            grid.pos[ip] = {
                wdata.radius * std::cos(iphi * 2.0 * M_PI / npts) wdata.radius *
                    std::sin(iphi * 2.0 * M_PI / npts),
                0.0};
            ++ip;
        }
    }
};

} // namespace ops
} // namespace ib
} // namespace amr_wind

#endif /* TAYLOR_GREEN_OPS_H */
