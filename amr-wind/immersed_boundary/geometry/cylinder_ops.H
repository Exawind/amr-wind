#ifndef SPHERE_OPS_H
#define SPHERE_OPS_H

#include "amr-wind/immersed_boundary/geometry/Cylinder.H"
#include "amr-wind/immersed_boundary/IBOps.H"
#include "amr-wind/immersed_boundary/IB.H"
#include "amr-wind/immersed_boundary/geometry/bluff_body_ops.H"

namespace amr_wind {
namespace ib {

namespace ops {

template <>
struct ReadInputsOp<Cylinder>
{
    void operator()(Cylinder::DataType& data, const utils::IBParser& pp)
    {
        auto& wdata = data.meta();
        auto& info = data.info();
        pp.get("num_points", wdata.num_pts);
        pp.get("center", wdata.center_loc);
        pp.get("radius", wdata.radius);
        pp.get("length", wdata.length);

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
struct InitDataOp<Cylinder>
{
    void operator()(Cylinder::DataType& data)
    {
        auto& wdata = data.meta();
        auto& grid = data.grid();
        int npts = wdata.num_pts;

        // Define the cylinder length
        amrex::Real dl = 2 * M_PI * wdata.radius / npts;
        // Compute the number of layers so that dl=dz
        int nlength = int(wdata.length / dl) + 1;
        dl = wdata.length / (nlength - 1);
        // Compute the number of radii so that dr=dl
        grid.resize(nlength * npts);
        amrex::Print()
            << "Total number of Lagrangian Markers for the Cylinder : "
            << nlength * npts << std::endl;

        wdata.dVl = dl * dl * dl;

        int ip = 0;
        for (int il = 0; il < nlength; ++il) {
            for (int iphi = 0; iphi < npts; ++iphi) {
                // clang-format off
                grid.pos[ip] = {wdata.center_loc[0] + wdata.radius* std::cos(iphi * 2.0 * M_PI / npts),
                                wdata.center_loc[1] - 0.5*wdata.length + il * dl,
                                wdata.center_loc[2] + wdata.radius* std::sin(iphi * 2.0 * M_PI / npts)};
                ++ip;
                // clang-formaton
            }
        }
        // This is no-slip only for a fixed sphere
        for (int i = 0; i < npts * nlength; ++i) {
            grid.forcing_vel[i] = {0.0, 0.0, 0.0};
        }
    
        auto& sim = data.sim();
        auto& mask_cell=sim.repo().get_int_field("mask_cell");
        auto& mask_node=sim.repo().get_int_field("mask_node");
        auto nlevels = sim.repo().num_active_levels();
        auto geom = sim.mesh().Geom();

        for (int lev = 0; lev < nlevels; ++lev) {
            const auto& problo = geom[lev].ProbLoArray();
            const auto& dx = geom[lev].CellSizeArray();
            amrex::Real mesh_size = std::cbrt(dx[0]*dx[1]*dx[2]);
            for (amrex::MFIter mfi(mask_cell(lev)); mfi.isValid(); ++mfi) {
                    const auto& vbx = mfi.growntilebox();
                    auto epsilon_cell = mask_cell(lev).array(mfi);
                amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                const amrex::Real x = problo[0] + (i + 0.5) * dx[0];
                const amrex::Real y = problo[1] + (j + 0.5) * dx[1];
                const amrex::Real z = problo[2] + (k + 0.5) * dx[2];
                
                amrex::Real r = std::sqrt((x-wdata.center_loc[0])*(x-wdata.center_loc[0])+
                                          (y-wdata.center_loc[1])*(y-wdata.center_loc[1])+
                                          (z-wdata.center_loc[2])*(z-wdata.center_loc[2]));
                if(r < wdata.radius + mesh_size){
                   epsilon_cell(i,j,k)=1;
                }

                });
            }

            for (amrex::MFIter mfi(mask_node(lev)); mfi.isValid(); ++mfi) {
                    const auto& vbx = mfi.growntilebox(1);
                    auto epsilon_node = mask_node(lev).array(mfi);
                amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                const amrex::Real x = problo[0] + i  * dx[0];
                const amrex::Real y = problo[1] + j  * dx[1];
                const amrex::Real z = problo[2] + k  * dx[2];
                
                amrex::Real r = std::sqrt((x-wdata.center_loc[0])*(x-wdata.center_loc[0])+
                                          (y-wdata.center_loc[1])*(y-wdata.center_loc[1])+
                                          (z-wdata.center_loc[2])*(z-wdata.center_loc[2]));
                if(r < wdata.radius + mesh_size){
                   epsilon_node(i,j,k)=1;
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
