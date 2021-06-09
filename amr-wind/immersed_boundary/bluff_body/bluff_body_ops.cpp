#include "amr-wind/immersed_boundary/bluff_body/bluff_body_ops.H"
#include "amr-wind/core/MultiParser.H"
#include "amr-wind/utilities/ncutils/nc_interface.H"
#include "amr-wind/utilities/io_utils.H"
// Used for mms
#include "amr-wind/physics/ConvectingTaylorVortex.H"

#include "AMReX_ParmParse.H"

namespace amr_wind {
namespace ib {
namespace bluff_body {

void read_inputs(
    BluffBodyBaseData& wdata, IBInfo&, const ::amr_wind::utils::MultiParser& pp)
{
    pp.query("has_wall_model", wdata.has_wall_model);
    pp.query("is_moving", wdata.is_moving);
    pp.query("is_mms", wdata.is_mms);
    pp.queryarr("vel_bc", wdata.vel_bc);
}

void init_data_structures(BluffBodyBaseData&) {}

void apply_mms_vel(CFDSim& sim)
{
    const int nlevels = sim.repo().num_active_levels();
    auto& mask_cell = sim.repo().get_int_field("mask_cell");
    auto& velocity = sim.repo().get_field("velocity");

    auto& m_conv_taylor_green =
        sim.physics_manager().get<ctv::ConvectingTaylorVortex>();

    const amrex::Real u0 = m_conv_taylor_green.get_u0();
    const amrex::Real v0 = m_conv_taylor_green.get_u0();
    const amrex::Real omega = m_conv_taylor_green.get_omega();

    amrex::Real t = sim.time().current_time();
    auto& geom = sim.mesh().Geom();

    for (int lev = 0; lev < nlevels; ++lev) {

        const auto& dx = geom[lev].CellSizeArray();
        const auto& problo = geom[lev].ProbLoArray();

        for (amrex::MFIter mfi(mask_cell(lev)); mfi.isValid(); ++mfi) {
            const auto& bx = mfi.growntilebox();
            auto epsilon_cell = mask_cell(lev).array(mfi);
            auto varr = velocity(lev).array(mfi);
            amrex::ParallelFor(
                bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    const amrex::Real x = problo[0] + (i + 0.5) * dx[0];
                    const amrex::Real y = problo[1] + (j + 0.5) * dx[1];
                    if (epsilon_cell(i, j, k) == 0) {
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

void apply_dirichlet_vel(CFDSim& sim, amrex::Vector<amrex::Real>& vel_bc)
{
    const int nlevels = sim.repo().num_active_levels();
    auto& mask_cell = sim.repo().get_int_field("mask_cell");
    auto& velocity = sim.repo().get_field("velocity");

    for (int lev = 0; lev < nlevels; ++lev) {
        for (amrex::MFIter mfi(mask_cell(lev)); mfi.isValid(); ++mfi) {
            const auto& bx = mfi.growntilebox();
            auto epsilon_cell = mask_cell(lev).array(mfi);
            auto varr = velocity(lev).array(mfi);
            amrex::ParallelFor(
                bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    if (epsilon_cell(i, j, k) == 0) {
                        varr(i, j, k, 0) = vel_bc[0];
                        varr(i, j, k, 1) = vel_bc[1];
                        varr(i, j, k, 2) = vel_bc[2];
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

} // namespace bluff_body
} // namespace ib
} // namespace amr_wind
