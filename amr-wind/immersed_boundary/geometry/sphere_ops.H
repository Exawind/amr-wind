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
        grid.resize(npts);

        int Ntheta, Nphi;
        // We discretise using spherical coordinates
        if (npts % 2 == 1) {
            Ntheta = (npts - 1) / 2;
        } else {
            Ntheta = npts / 2;
        }
        Nphi = npts - Ntheta;
        int ip = 0;
        for (int itheta = 0; itheta < Ntheta; ++itheta) {
            for (int iphi = 0; iphi < Nphi; ++iphi) {
                grid.pos[ip] = {
                    wdata.radius * std::sin(iphi * 2.0 * M_PI / Nphi) *
                        std::cos(itheta * 2.0 * M_PI / Ntheta),
                    wdata.radius * std::sin(iphi * 2.0 * M_PI / Nphi) *
                        std::sin(itheta * 2.0 * M_PI / Ntheta),
                    wdata.radius * std::cos(iphi * 2.0 * M_PI / Nphi)};
            }
            ++ip;
        }
    }
};

} // namespace ops
} // namespace ib
} // namespace amr_wind

#endif /* SPHERE_OPS_H */
