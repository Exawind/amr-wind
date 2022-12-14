#include "amr-wind/immersed_boundary/bluff_body/bluff_body_ops.H"
#include "amr-wind/core/MultiParser.H"
#include "amr-wind/utilities/ncutils/nc_interface.H"
#include "amr-wind/utilities/io_utils.H"

#include "amr-wind/fvm/gradient.H"
#include "amr-wind/core/field_ops.H"

// Used for mms
#include "amr-wind/physics/ConvectingTaylorVortex.H"

#include "AMReX_ParmParse.H"

namespace amr_wind::ib::bluff_body {

void read_inputs(
    BluffBodyBaseData& wdata,
    IBInfo& /*unused*/,
    const ::amr_wind::utils::MultiParser& pp)
{
    pp.query("has_wall_model", wdata.has_wall_model);
    pp.query("is_moving", wdata.is_moving);
    pp.query("is_mms", wdata.is_mms);
    pp.queryarr("vel_bc", wdata.vel_bc);
}

void init_data_structures(BluffBodyBaseData& /*unused*/) {}

void apply_mms_vel(CFDSim& sim)
{
    const int nlevels = sim.repo().num_active_levels();

    const auto& levelset = sim.repo().get_field("ib_levelset");
    // cppcheck-suppress constVariable
    auto& velocity = sim.repo().get_field("velocity");
    auto& m_conv_taylor_green =
        sim.physics_manager().get<ctv::ConvectingTaylorVortex>();

    const amrex::Real u0 = m_conv_taylor_green.get_u0();
    const amrex::Real v0 = m_conv_taylor_green.get_v0();
    const amrex::Real omega = m_conv_taylor_green.get_omega();

    amrex::Real t = sim.time().new_time();
    auto& geom = sim.mesh().Geom();

    for (int lev = 0; lev < nlevels; ++lev) {

        const auto& dx = geom[lev].CellSizeArray();
        const auto& problo = geom[lev].ProbLoArray();

        for (amrex::MFIter mfi(levelset(lev)); mfi.isValid(); ++mfi) {
            const auto& bx = mfi.growntilebox();
            auto phi = levelset(lev).array(mfi);
            auto varr = velocity(lev).array(mfi);
            amrex::ParallelFor(
                bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    const amrex::Real x = problo[0] + (i + 0.5) * dx[0];
                    const amrex::Real y = problo[1] + (j + 0.5) * dx[1];

                    if (phi(i, j, k) <= 0) {
                        varr(i, j, k, 0) =
                            u0 - std::cos(utils::pi() * (x - u0 * t)) *
                                     std::sin(utils::pi() * (y - v0 * t)) *
                                     std::exp(-2.0 * omega * t);
                        varr(i, j, k, 1) =
                            v0 + std::sin(utils::pi() * (x - u0 * t)) *
                                     std::cos(utils::pi() * (y - v0 * t)) *
                                     std::exp(-2.0 * omega * t);
                        varr(i, j, k, 2) = 0.0;
                    }
                });
        }
    }
}

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
        // Defining the "ghost-cell" band distance
        amrex::Real phi_b = std::cbrt(dx[0] * dx[1] * dx[2]);

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
                    // Pure solid-body points
                    if (phi_arr(i, j, k) < -phi_b) {
                        varr(i, j, k, 0) = velx;
                        varr(i, j, k, 1) = vely;
                        varr(i, j, k, 2) = velz;
                        norm_arr(i, j, k, 0) = 0.;
                        norm_arr(i, j, k, 1) = 0.;
                        norm_arr(i, j, k, 2) = 0.;
                        // This determines the ghost-cells
                    } else if (
                        phi_arr(i, j, k) < 0 && phi_arr(i, j, k) >= -phi_b) {
                        // For this particular ghost-cell find the
                        // body-intercept (BI) point and image-point (IP)
                        // First define the ghost cell point
                        // amrex::Real x_GC = problo[0] + (i + 0.5) * dx[0];
                        // amrex::Real y_GC = problo[1] + (j + 0.5) * dx[1];
                        // amrex::Real z_GC = problo[2] + (k + 0.5) * dx[2];
                        // Find the "image-points"
                        varr(i, j, k, 0) = velx;
                        varr(i, j, k, 1) = vely;
                        varr(i, j, k, 2) = velz;
                    } else {
                        norm_arr(i, j, k, 0) = 0.;
                        norm_arr(i, j, k, 1) = 0.;
                        norm_arr(i, j, k, 2) = 0.;
                    }
                });
        }
    }
}

void prepare_netcdf_file(
    const std::string& ncfile,
    const BluffBodyBaseData& meta,
    const IBInfo& info)
{
    amrex::ignore_unused(ncfile, meta, info);
}

void write_netcdf(
    const std::string& ncfile,
    const BluffBodyBaseData& meta,
    const IBInfo& info,
    const amrex::Real time)
{
    amrex::ignore_unused(ncfile, meta, info, time);
}

} // namespace amr_wind::ib::bluff_body
