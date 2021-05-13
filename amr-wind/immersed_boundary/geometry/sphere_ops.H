#ifndef SPHERE_OPS_H
#define SPHERE_OPS_H

#include "amr-wind/immersed_boundary/geometry/Sphere.H"
#include "amr-wind/immersed_boundary/IBOps.H"
#include "amr-wind/immersed_boundary/geometry/bluff_body_ops.H"

namespace amr_wind {
namespace ib {

namespace ops {

template <>
struct ReadInputsOp<Sphere>
{
    void operator()(Sphere::DataType& data, const utils::IBParser& pp)
    {
        auto& wdata = data.meta();
        auto& info = data.info();
        pp.get("num_points", wdata.num_pts);
        pp.get("center", wdata.center_loc);
        pp.get("radius", wdata.radius);
        amrex::Real search_radius = 1.5 * wdata.radius;
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
        auto& grid = data.grid();
        int npts = wdata.num_pts;
        grid.resize(npts * npts);

        int ip = 0;
        for (int itheta = 0; itheta < npts; ++itheta) {
            for (int iphi = 0; iphi < npts; ++iphi) {
                grid.pos[ip] = {
                    wdata.radius * std::sin(iphi * 2.0 * M_PI / npts) *
                        std::cos(itheta * 2.0 * M_PI / npts),
                    wdata.radius * std::sin(iphi * 2.0 * M_PI / npts) *
                        std::sin(itheta * 2.0 * M_PI / npts),
                    wdata.radius * std::cos(iphi * 2.0 * M_PI / npts)};
                ++ip;
            }
        }
        // This is no-slip only for a fixed sphere
        for (int i = 0; i < npts * npts; ++i) {
            grid.forcing_vel[i] = {0.0, 0.0, 0.0};
        }
    }
};

} // namespace ops
} // namespace ib
} // namespace amr_wind

#endif /* SPHERE_OPS_H */
