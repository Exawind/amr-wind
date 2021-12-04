#include <cmath>

#include "amr-wind/mesh_mapping_models/ChannelFlowScaling.H"
#include "amr-wind/CFDSim.H"

#include "AMReX_ParmParse.H"

namespace amr_wind {

ChannelFlowScaling::ChannelFlowScaling(const CFDSim& sim)
    : m_mesh_scale_fac_cc(sim.repo().get_field("mesh_scaling_factor_cc"))
    , m_mesh_scale_fac_nd(sim.repo().get_field("mesh_scaling_factor_nd"))
{
    {
        amrex::ParmParse pp("ChannelFlowScaling");
        pp.queryarr("beta", m_beta, 0, AMREX_SPACEDIM);
        pp.queryarr("do_map", m_map, 0, AMREX_SPACEDIM);
    }

    {
        amrex::ParmParse pp("amr");
        pp.getarr("n_cell", n_cell, 0, AMREX_SPACEDIM);
    }
}

/** Construct the mesh mapping field
 */
void ChannelFlowScaling::create_map(int lev, const amrex::Geometry& geom)
{
    const auto beta = m_beta;
    const auto do_map = m_map;
    const auto num_cell = n_cell;

    const auto& dx = geom.CellSizeArray();
    const auto& prob_lo = geom.ProbLoArray();
    const auto& prob_hi = geom.ProbHiArray();

    amrex::Vector<amrex::Real> len_cc{
        {(prob_hi[0] - dx[0] / 2) - (prob_lo[0] + dx[0] / 2),
         (prob_hi[1] - dx[1] / 2) - (prob_lo[1] + dx[1] / 2),
         (prob_hi[2] - dx[2] / 2) - (prob_lo[2] + dx[2] / 2)}};
    amrex::Vector<amrex::Real> len_nd{
        {prob_hi[0] - prob_lo[0], prob_hi[1] - prob_lo[1],
         prob_hi[2] - prob_lo[2]}};

    for (amrex::MFIter mfi(m_mesh_scale_fac_cc(lev)); mfi.isValid(); ++mfi) {

        const auto& bx = mfi.growntilebox();
        amrex::Array4<amrex::Real> const& scale_fac_cc =
            m_mesh_scale_fac_cc(lev).array(mfi);
        amrex::ParallelFor(
            bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];
                amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1];
                amrex::Real z = prob_lo[2] + (k + 0.5) * dx[2];

                amrex::Real fac_x =
                    beta[0] *
                    (1 -
                     std::pow(
                         std::tanh(
                             beta[0] * (1 - 2 * (x - prob_lo[0]) / len_cc[0])),
                         2)) /
                    std::tanh(beta[0]);
                amrex::Real fac_y =
                    beta[1] *
                    (1 -
                     std::pow(
                         std::tanh(
                             beta[1] * (1 - 2 * (y - prob_lo[1]) / len_cc[1])),
                         2)) /
                    std::tanh(beta[1]);
                amrex::Real fac_z =
                    beta[2] *
                    (1 -
                     std::pow(
                         std::tanh(
                             beta[2] * (1 - 2 * (z - prob_lo[2]) / len_cc[2])),
                         2)) /
                    std::tanh(beta[2]);

                scale_fac_cc(i, j, k, 0) =
                    (do_map[0] && (i >= 0) && (i < num_cell[0])) ? fac_x : 1.0;
                scale_fac_cc(i, j, k, 1) =
                    (do_map[1] && (j >= 0) && (j < num_cell[1])) ? fac_y : 1.0;
                scale_fac_cc(i, j, k, 2) =
                    (do_map[2] && (k >= 0) && (k < num_cell[2])) ? fac_z : 1.0;
            });

        const auto& nbx = mfi.grownnodaltilebox();
        amrex::Array4<amrex::Real> const& scale_fac_nd =
            m_mesh_scale_fac_nd(lev).array(mfi);
        amrex::ParallelFor(
            nbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                amrex::Real x = prob_lo[0] + i * dx[0];
                amrex::Real y = prob_lo[1] + j * dx[1];
                amrex::Real z = prob_lo[2] + k * dx[2];

                amrex::Real fac_x =
                    beta[0] *
                    (1 -
                     std::pow(
                         std::tanh(
                             beta[0] * (1 - 2 * (x - prob_lo[0]) / len_nd[0])),
                         2)) /
                    std::tanh(beta[0]);
                amrex::Real fac_y =
                    beta[1] *
                    (1 -
                     std::pow(
                         std::tanh(
                             beta[1] * (1 - 2 * (y - prob_lo[1]) / len_nd[1])),
                         2)) /
                    std::tanh(beta[1]);
                amrex::Real fac_z =
                    beta[2] *
                    (1 -
                     std::pow(
                         std::tanh(
                             beta[2] * (1 - 2 * (z - prob_lo[2]) / len_nd[2])),
                         2)) /
                    std::tanh(beta[2]);

                scale_fac_nd(i, j, k, 0) =
                    (do_map[0] && (i >= 0) && (i <= num_cell[0])) ? fac_x : 1.0;
                scale_fac_nd(i, j, k, 1) =
                    (do_map[1] && (j >= 0) && (j <= num_cell[1])) ? fac_y : 1.0;
                scale_fac_nd(i, j, k, 2) =
                    (do_map[2] && (k >= 0) && (k <= num_cell[2])) ? fac_z : 1.0;
            });
    }

    // TODO: Call fill patch operators
}

} // namespace amr_wind
