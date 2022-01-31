#include "amr-wind/immersed_boundary/complex_terrain/complex_terrain_ops.H"
#include "amr-wind/core/MultiParser.H"
#include "amr-wind/utilities/ncutils/nc_interface.H"
#include "amr-wind/utilities/io_utils.H"

#include "amr-wind/fvm/gradient.H"
#include "amr-wind/core/field_ops.H"

// Used for mms
#include "amr-wind/physics/ConvectingTaylorVortex.H"

#include "AMReX_ParmParse.H"

namespace amr_wind {
namespace ib {
namespace complex_terrain {

void read_inputs(
    ComplexTerrainBaseData& wdata,
    IBInfo& /*unused*/,
    const ::amr_wind::utils::MultiParser& pp)
{
    pp.query("has_wall_model", wdata.has_wall_model);
}

void init_data_structures(ComplexTerrainBaseData& /*unused*/) {}

void apply_dirichlet_vel(CFDSim& sim, const amrex::Vector<amrex::Real>& vel_bc)
{
    const int nlevels = sim.repo().num_active_levels();
    auto& geom = sim.mesh().Geom();
    // cppcheck-suppress constVariable
    auto& velocity = sim.repo().get_field("velocity");
    auto& levelset = sim.repo().get_field("ib_levelset");
    levelset.fillpatch(sim.time().current_time());
    auto& normal = sim.repo().get_field("ib_normal");
    fvm::gradient(normal, levelset);
    field_ops::normalize(normal);
    normal.fillpatch(sim.time().current_time());

    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& dx = geom[lev].CellSizeArray();
        // const auto& problo = geom[lev].ProbLoArray();
        // Defining the band distance
        amrex::Real phi_b = 2. * std::cbrt(dx[0] * dx[1] * dx[2]);

        for (amrex::MFIter mfi(levelset(lev)); mfi.isValid(); ++mfi) {
            const auto& bx = mfi.tilebox();
            auto varr = velocity(lev).array(mfi);
            auto phi_arr = levelset(lev).array(mfi);
            auto norm_arr = normal(lev).array(mfi);

            amrex::Real velx = vel_bc[0];
            amrex::Real vely = vel_bc[1];
            amrex::Real velz = vel_bc[2];

            amrex::ParallelFor(
                bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    if (phi_arr(i, j, k) < 0) {
                        varr(i, j, k, 0) = velx;
                        varr(i, j, k, 1) = vely;
                        varr(i, j, k, 2) = velz;
                        norm_arr(i, j, k, 0) = 0.;
                        norm_arr(i, j, k, 1) = 0.;
                        norm_arr(i, j, k, 2) = 0.;
                    } else if (phi_arr(i, j, k) > phi_b) {
                        norm_arr(i, j, k, 0) = 0.;
                        norm_arr(i, j, k, 1) = 0.;
                        norm_arr(i, j, k, 2) = 0.;
                    }
                });
        }
    }
}

void compute_wall_stresses(
    CFDSim& sim, const amrex::Vector<amrex::Real>& tau_wall)
{}

void prepare_netcdf_file(
    const std::string& ncfile,
    const ComplexTerrainBaseData& meta,
    const IBInfo& info)
{
    amrex::ignore_unused(ncfile, meta, info);
}

void write_netcdf(
    const std::string& ncfile,
    const ComplexTerrainBaseData& meta,
    const IBInfo& info,
    const amrex::Real time)
{
    amrex::ignore_unused(ncfile, meta, info, time);
}

} // namespace complex_terrain
} // namespace ib
} // namespace amr_wind
