#ifndef RECTANGULAR_PRISM_OPS_H
#define RECTANGULAR_PRISM_OPS_H

#include "amr-wind/immersed_boundary/geometry/RectangularPrism.H"
#include "amr-wind/immersed_boundary/IBOps.H"
#include "amr-wind/immersed_boundary/IB.H"
#include "amr-wind/immersed_boundary/geometry/bluff_body_ops.H"

namespace amr_wind {
namespace ib {

namespace ops {

template <>
struct ReadInputsOp<RectangularPrism>
{
    void operator()(RectangularPrism::DataType& data, const utils::IBParser& pp)
    {
        auto& wdata = data.meta();
        auto& info = data.info();
        pp.get("num_points", wdata.num_pts);
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
struct InitDataOp<RectangularPrism>
{
    void operator()(RectangularPrism::DataType& data)
    {
        auto& wdata = data.meta();
        auto& grid = data.grid();
        int npts = wdata.num_pts;

        // Define the prisgm length
        amrex::Real dl = wdata.length / npts;
        int nlength = int(wdata.length / dl) + 1;

        amrex::Real dh = wdata.height / npts;
        int nheight = int(wdata.height / dh) + 1;

        amrex::Real dw = wdata.width / npts;
        int nwidth = int(wdata.width / dw) + 1;

        // Compute the number of radii so that dr=dl
        grid.resize(npts * npts * npts);
        amrex::Print() << "Total number of Lagrangian Markers for the "
                       << RectangularPrism::identifier() << " : "
                       << npts * npts * npts << std::endl;

        wdata.dVl = dl * dh * dw;

        int ip = 0;
        for (int il = 0; il < nlength; ++il) {
            for (int iw = 0; iw < nwidth; ++iw) {
                for (int ih = 0; ih < nheight; ++ih) {

                    // clang-format off
                grid.pos[ip] = {wdata.center_loc[0] - 0.5*wdata.length + il * dl,
                                wdata.center_loc[1] - 0.5*wdata.width  + iw * dw,
                                wdata.center_loc[2] - 0.5*wdata.height + ih * dh};
                ++ip;
                // clang-formaton
                }
            }
        }
        // This is no-slip only for a fixed sphere
        for (int i = 0; i < npts * nlength; ++i) {
            grid.forcing_vel[i] = {0.0, 0.0, 0.0};
        }
        
        auto& sim = data.sim();
        auto nlevels = sim.repo().num_active_levels();
        auto geom = sim.mesh().Geom();
        for (int lev = 0; lev < nlevels; ++lev) {
            const auto& dx = geom[lev].CellSizeArray();
            amrex::Print() << "Lagrangian marker volume to mesh volume ratio at level "
                           << lev << " : "<<wdata.dVl/(dx[0]*dx[1]*dx[2])<< std::endl;

        }

    }
};

} // namespace ops
} // namespace ib
} // namespace amr_wind

#endif /* RectangularPrism_OPS_H */
