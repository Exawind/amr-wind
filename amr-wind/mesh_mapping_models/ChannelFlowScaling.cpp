#include <cmath>

#include "amr-wind/mesh_mapping_models/ChannelFlowScaling.H"
#include "amr-wind/CFDSim.H"

#include "AMReX_ParmParse.H"

namespace amr_wind {

ChannelFlowScaling::ChannelFlowScaling(const CFDSim& sim)
    : m_mesh_scale_fac_cc(sim.repo().get_field("mesh_scaling_factor_cc"))
    , m_mesh_scale_fac_nd(sim.repo().get_field("mesh_scaling_factor_nd"))
    , m_mesh_scale_fac_xf(sim.repo().get_field("mesh_scaling_factor_xf"))
    , m_mesh_scale_fac_yf(sim.repo().get_field("mesh_scaling_factor_yf"))
    , m_mesh_scale_fac_zf(sim.repo().get_field("mesh_scaling_factor_zf"))
    , m_non_uniform_coord_cc(sim.repo().get_field("non_uniform_coord_cc"))
    , m_non_uniform_coord_nd(sim.repo().get_field("non_uniform_coord_nd"))
{
    amrex::ParmParse pp("ChannelFlowScaling");
    pp.queryarr("beta", m_beta, 0, AMREX_SPACEDIM);
    pp.queryarr("do_map", m_map, 0, AMREX_SPACEDIM);
}

/** Construct the mesh mapping field
 */
void ChannelFlowScaling::create_map(int lev, const amrex::Geometry& geom)
{
    create_cell_node_map(lev, geom);
    create_face_map(lev, geom);
    create_non_uniform_mesh(lev, geom);
}

/** Construct the mesh mapping field on cell centers and nodes
 */
void ChannelFlowScaling::create_cell_node_map(
    int lev, const amrex::Geometry& geom)
{
    const auto beta = m_beta;
    const auto do_map = m_map;
    const auto eps = m_eps;

    const auto& dx = geom.CellSizeArray();
    const auto& prob_lo = geom.ProbLoArray();
    const auto& prob_hi = geom.ProbHiArray();

    amrex::Vector<amrex::Real> len{
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
                    (1 - std::pow(
                             std::tanh(
                                 beta[0] * (1 - 2 * (x - prob_lo[0]) / len[0])),
                             2)) /
                    std::tanh(beta[0]);
                amrex::Real fac_y =
                    beta[1] *
                    (1 - std::pow(
                             std::tanh(
                                 beta[1] * (1 - 2 * (y - prob_lo[1]) / len[1])),
                             2)) /
                    std::tanh(beta[1]);
                amrex::Real fac_z =
                    beta[2] *
                    (1 - std::pow(
                             std::tanh(
                                 beta[2] * (1 - 2 * (z - prob_lo[2]) / len[2])),
                             2)) /
                    std::tanh(beta[2]);

                bool in_domain =
                    ((x > prob_lo[0]) && (x < prob_hi[0]) && (y > prob_lo[1]) &&
                     (y < prob_hi[1]) && (z > prob_lo[2]) && (z < prob_hi[2]));

                scale_fac_cc(i, j, k, 0) =
                    (do_map[0] && in_domain) ? fac_x : 1.0;
                scale_fac_cc(i, j, k, 1) =
                    (do_map[1] && in_domain) ? fac_y : 1.0;
                scale_fac_cc(i, j, k, 2) =
                    (do_map[2] && in_domain) ? fac_z : 1.0;
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
                    (1 - std::pow(
                             std::tanh(
                                 beta[0] * (1 - 2 * (x - prob_lo[0]) / len[0])),
                             2)) /
                    std::tanh(beta[0]);
                amrex::Real fac_y =
                    beta[1] *
                    (1 - std::pow(
                             std::tanh(
                                 beta[1] * (1 - 2 * (y - prob_lo[1]) / len[1])),
                             2)) /
                    std::tanh(beta[1]);
                amrex::Real fac_z =
                    beta[2] *
                    (1 - std::pow(
                             std::tanh(
                                 beta[2] * (1 - 2 * (z - prob_lo[2]) / len[2])),
                             2)) /
                    std::tanh(beta[2]);

                bool in_domain =
                    ((x >= prob_lo[0] - eps) && (x <= prob_hi[0] + eps) &&
                     (y >= prob_lo[1] - eps) && (y <= prob_hi[1] + eps) &&
                     (z >= prob_lo[2] - eps) && (z <= prob_hi[2] + eps));

                scale_fac_nd(i, j, k, 0) =
                    (do_map[0] && in_domain) ? fac_x : 1.0;
                scale_fac_nd(i, j, k, 1) =
                    (do_map[1] && in_domain) ? fac_y : 1.0;
                scale_fac_nd(i, j, k, 2) =
                    (do_map[2] && in_domain) ? fac_z : 1.0;
            });
    }

    // TODO: Call fill patch operators
}

/** Construct the mesh mapping field on cell faces
 */
void ChannelFlowScaling::create_face_map(int lev, const amrex::Geometry& geom)
{
    const auto beta = m_beta;
    const auto do_map = m_map;
    const auto eps = m_eps;

    const auto& dx = geom.CellSizeArray();
    const auto& prob_lo = geom.ProbLoArray();
    const auto& prob_hi = geom.ProbHiArray();

    amrex::Vector<amrex::Real> len{
        {prob_hi[0] - prob_lo[0], prob_hi[1] - prob_lo[1],
         prob_hi[2] - prob_lo[2]}};

    for (amrex::MFIter mfi(m_mesh_scale_fac_xf(lev)); mfi.isValid(); ++mfi) {

        const auto& bx = mfi.growntilebox();
        amrex::Array4<amrex::Real> const& scale_fac_xf =
            m_mesh_scale_fac_xf(lev).array(mfi);
        amrex::ParallelFor(
            bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                amrex::Real x = prob_lo[0] + i * dx[0];
                amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1];
                amrex::Real z = prob_lo[2] + (k + 0.5) * dx[2];

                amrex::Real fac_x =
                    beta[0] *
                    (1 - std::pow(
                             std::tanh(
                                 beta[0] * (1 - 2 * (x - prob_lo[0]) / len[0])),
                             2)) /
                    std::tanh(beta[0]);
                amrex::Real fac_y =
                    beta[1] *
                    (1 - std::pow(
                             std::tanh(
                                 beta[1] * (1 - 2 * (y - prob_lo[1]) / len[1])),
                             2)) /
                    std::tanh(beta[1]);
                amrex::Real fac_z =
                    beta[2] *
                    (1 - std::pow(
                             std::tanh(
                                 beta[2] * (1 - 2 * (z - prob_lo[2]) / len[2])),
                             2)) /
                    std::tanh(beta[2]);

                bool in_domain =
                    ((x >= prob_lo[0] - eps) && (x <= prob_hi[0] + eps) &&
                     (y > prob_lo[1]) && (y < prob_hi[1]) && (z > prob_lo[2]) &&
                     (z < prob_hi[2]));

                scale_fac_xf(i, j, k, 0) =
                    (do_map[0] && in_domain) ? fac_x : 1.0;
                scale_fac_xf(i, j, k, 1) =
                    (do_map[1] && in_domain) ? fac_y : 1.0;
                scale_fac_xf(i, j, k, 2) =
                    (do_map[2] && in_domain) ? fac_z : 1.0;
            });
    }

    for (amrex::MFIter mfi(m_mesh_scale_fac_yf(lev)); mfi.isValid(); ++mfi) {

        const auto& bx = mfi.growntilebox();
        amrex::Array4<amrex::Real> const& scale_fac_yf =
            m_mesh_scale_fac_yf(lev).array(mfi);
        amrex::ParallelFor(
            bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];
                amrex::Real y = prob_lo[1] + j * dx[1];
                amrex::Real z = prob_lo[2] + (k + 0.5) * dx[2];

                amrex::Real fac_x =
                    beta[0] *
                    (1 - std::pow(
                             std::tanh(
                                 beta[0] * (1 - 2 * (x - prob_lo[0]) / len[0])),
                             2)) /
                    std::tanh(beta[0]);
                amrex::Real fac_y =
                    beta[1] *
                    (1 - std::pow(
                             std::tanh(
                                 beta[1] * (1 - 2 * (y - prob_lo[1]) / len[1])),
                             2)) /
                    std::tanh(beta[1]);
                amrex::Real fac_z =
                    beta[2] *
                    (1 - std::pow(
                             std::tanh(
                                 beta[2] * (1 - 2 * (z - prob_lo[2]) / len[2])),
                             2)) /
                    std::tanh(beta[2]);

                bool in_domain =
                    ((x > prob_lo[0]) && (x < prob_hi[0]) &&
                     (y >= prob_lo[1] - eps) && (y <= prob_hi[1] + eps) &&
                     (z > prob_lo[2]) && (z < prob_hi[2]));

                scale_fac_yf(i, j, k, 0) =
                    (do_map[0] && in_domain) ? fac_x : 1.0;
                scale_fac_yf(i, j, k, 1) =
                    (do_map[1] && in_domain) ? fac_y : 1.0;
                scale_fac_yf(i, j, k, 2) =
                    (do_map[2] && in_domain) ? fac_z : 1.0;
            });
    }

    for (amrex::MFIter mfi(m_mesh_scale_fac_zf(lev)); mfi.isValid(); ++mfi) {

        const auto& bx = mfi.growntilebox();
        amrex::Array4<amrex::Real> const& scale_fac_zf =
            m_mesh_scale_fac_zf(lev).array(mfi);
        amrex::ParallelFor(
            bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];
                amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1];
                amrex::Real z = prob_lo[2] + k * dx[2];

                amrex::Real fac_x =
                    beta[0] *
                    (1 - std::pow(
                             std::tanh(
                                 beta[0] * (1 - 2 * (x - prob_lo[0]) / len[0])),
                             2)) /
                    std::tanh(beta[0]);
                amrex::Real fac_y =
                    beta[1] *
                    (1 - std::pow(
                             std::tanh(
                                 beta[1] * (1 - 2 * (y - prob_lo[1]) / len[1])),
                             2)) /
                    std::tanh(beta[1]);
                amrex::Real fac_z =
                    beta[2] *
                    (1 - std::pow(
                             std::tanh(
                                 beta[2] * (1 - 2 * (z - prob_lo[2]) / len[2])),
                             2)) /
                    std::tanh(beta[2]);

                bool in_domain =
                    ((x > prob_lo[0]) && (x < prob_hi[0]) && (y > prob_lo[1]) &&
                     (y < prob_hi[1]) && (z >= prob_lo[2] - eps) &&
                     (z <= prob_hi[2] + eps));

                scale_fac_zf(i, j, k, 0) =
                    (do_map[0] && in_domain) ? fac_x : 1.0;
                scale_fac_zf(i, j, k, 1) =
                    (do_map[1] && in_domain) ? fac_y : 1.0;
                scale_fac_zf(i, j, k, 2) =
                    (do_map[2] && in_domain) ? fac_z : 1.0;
            });
    }

    // TODO: Call fill patch operators
}

/** Construct the non-uniform mesh field
 */
void ChannelFlowScaling::create_non_uniform_mesh(
    int lev, const amrex::Geometry& geom)
{
    amrex::Vector<amrex::Real> probhi_physical{{0.0, 0.0, 0.0}};
    {
        amrex::ParmParse pp("geometry");
        if (pp.contains("prob_hi_physical")) {
            pp.getarr("prob_hi_physical", probhi_physical);
        } else {
            for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                probhi_physical[d] = geom.ProbHiArray()[d];
            }
        }
    }

    const auto beta = m_beta;
    const auto do_map = m_map;
    const auto eps = m_eps;

    const auto& dx = geom.CellSizeArray();
    const auto& prob_lo = geom.ProbLoArray();
    const auto& prob_hi = geom.ProbHiArray();

    amrex::Vector<amrex::Real> len{
        {probhi_physical[0] - prob_lo[0], probhi_physical[1] - prob_lo[1],
         probhi_physical[2] - prob_lo[2]}};

    for (amrex::MFIter mfi(m_non_uniform_coord_cc(lev)); mfi.isValid(); ++mfi) {

        const auto& bx = mfi.growntilebox();
        amrex::Array4<amrex::Real> const& nu_coord_cc =
            m_non_uniform_coord_cc(lev).array(mfi);
        amrex::ParallelFor(
            bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];
                amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1];
                amrex::Real z = prob_lo[2] + (k + 0.5) * dx[2];

                amrex::Real x_non_uni =
                    prob_lo[0] +
                    len[0] / 2 *
                        (1 -
                         std::tanh(
                             beta[0] * (1 - 2 * (x - prob_lo[0]) / len[0])) /
                             std::tanh(beta[0]));
                amrex::Real y_non_uni =
                    prob_lo[1] +
                    len[1] / 2 *
                        (1 -
                         std::tanh(
                             beta[1] * (1 - 2 * (y - prob_lo[1]) / len[1])) /
                             std::tanh(beta[1]));
                amrex::Real z_non_uni =
                    prob_lo[2] +
                    len[2] / 2 *
                        (1 -
                         std::tanh(
                             beta[2] * (1 - 2 * (z - prob_lo[2]) / len[2])) /
                             std::tanh(beta[2]));

                bool in_domain =
                    ((x > prob_lo[0]) && (x < prob_hi[0]) && (y > prob_lo[1]) &&
                     (y < prob_hi[1]) && (z > prob_lo[2]) && (z < prob_hi[2]));

                nu_coord_cc(i, j, k, 0) =
                    (do_map[0] && in_domain) ? x_non_uni : x;
                nu_coord_cc(i, j, k, 1) =
                    (do_map[1] && in_domain) ? y_non_uni : y;
                nu_coord_cc(i, j, k, 2) =
                    (do_map[2] && in_domain) ? z_non_uni : z;
            });

        const auto& nbx = mfi.grownnodaltilebox();
        amrex::Array4<amrex::Real> const& nu_coord_nd =
            m_non_uniform_coord_nd(lev).array(mfi);
        amrex::ParallelFor(
            nbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                amrex::Real x = prob_lo[0] + i * dx[0];
                amrex::Real y = prob_lo[1] + j * dx[1];
                amrex::Real z = prob_lo[2] + k * dx[2];

                amrex::Real x_non_uni =
                    prob_lo[0] +
                    len[0] / 2 *
                        (1 -
                         std::tanh(
                             beta[0] * (1 - 2 * (x - prob_lo[0]) / len[0])) /
                             std::tanh(beta[0]));
                amrex::Real y_non_uni =
                    prob_lo[1] +
                    len[1] / 2 *
                        (1 -
                         std::tanh(
                             beta[1] * (1 - 2 * (y - prob_lo[1]) / len[1])) /
                             std::tanh(beta[1]));
                amrex::Real z_non_uni =
                    prob_lo[2] +
                    len[2] / 2 *
                        (1 -
                         std::tanh(
                             beta[2] * (1 - 2 * (z - prob_lo[2]) / len[2])) /
                             std::tanh(beta[2]));

                bool in_domain =
                    ((x >= prob_lo[0] - eps) && (x <= prob_hi[0] + eps) &&
                     (y >= prob_lo[1] - eps) && (y <= prob_hi[1] + eps) &&
                     (z >= prob_lo[2] - eps) && (z <= prob_hi[2] + eps));

                nu_coord_nd(i, j, k, 0) =
                    (do_map[0] && in_domain) ? x_non_uni : x;
                nu_coord_nd(i, j, k, 1) =
                    (do_map[1] && in_domain) ? y_non_uni : y;
                nu_coord_nd(i, j, k, 2) =
                    (do_map[2] && in_domain) ? z_non_uni : z;
            });
    }
}

} // namespace amr_wind
