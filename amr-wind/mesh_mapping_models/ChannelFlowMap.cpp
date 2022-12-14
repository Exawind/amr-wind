#include <cmath>

#include "amr-wind/mesh_mapping_models/ChannelFlowMap.H"

#include "AMReX_ParmParse.H"

namespace amr_wind::channel_map {

namespace {

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real eval_fac(
    const amrex::Real x,
    const amrex::Real beta,
    const amrex::Real prob_lo,
    const amrex::Real len)
{
    return (beta == 0.0)
               ? 1.0
               : (beta *
                  (1 -
                   std::pow(
                       std::tanh(beta * (1 - 2 * (x - prob_lo) / len)), 2)) /
                  std::tanh(beta));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real eval_coord(
    const amrex::Real x,
    const amrex::Real beta,
    const amrex::Real prob_lo,
    const amrex::Real len)
{
    return (beta == 0.0)
               ? x
               : (prob_lo +
                  len / 2 *
                      (1 - std::tanh(beta * (1 - 2 * (x - prob_lo) / len)) /
                               std::tanh(beta)));
}

} // namespace

ChannelFlowMap::ChannelFlowMap()
{
    amrex::ParmParse pp("ChannelFlowMap");
    pp.queryarr("beta", m_beta, 0, AMREX_SPACEDIM);
}

/** Construct the mesh mapping field
 */
void ChannelFlowMap::create_map(int lev, const amrex::Geometry& geom)
{
    create_cell_node_map(lev, geom);
    create_face_map(lev, geom);
    create_non_uniform_mesh(lev, geom);
}

/** Construct the mesh mapping field on cell centers and nodes
 */
void ChannelFlowMap::create_cell_node_map(int lev, const amrex::Geometry& geom)
{
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> beta{
        {m_beta[0], m_beta[1], m_beta[2]}};
    const auto eps = m_eps;

    const auto& dx = geom.CellSizeArray();
    const auto& prob_lo = geom.ProbLoArray();
    const auto& prob_hi = geom.ProbHiArray();

    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> len{
        {prob_hi[0] - prob_lo[0], prob_hi[1] - prob_lo[1],
         prob_hi[2] - prob_lo[2]}};

    for (amrex::MFIter mfi((*m_mesh_scale_fac_cc)(lev)); mfi.isValid(); ++mfi) {

        const auto& bx = mfi.growntilebox();
        amrex::Array4<amrex::Real> const& scale_fac_cc =
            (*m_mesh_scale_fac_cc)(lev).array(mfi);
        amrex::Array4<amrex::Real> const& scale_detJ_cc =
            (*m_mesh_scale_detJ_cc)(lev).array(mfi);
        amrex::ParallelFor(
            bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];
                amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1];
                amrex::Real z = prob_lo[2] + (k + 0.5) * dx[2];

                amrex::Real fac_x = eval_fac(x, beta[0], prob_lo[0], len[0]);
                amrex::Real fac_y = eval_fac(y, beta[1], prob_lo[1], len[1]);
                amrex::Real fac_z = eval_fac(z, beta[2], prob_lo[2], len[2]);

                bool in_domain =
                    ((x > prob_lo[0]) && (x < prob_hi[0]) && (y > prob_lo[1]) &&
                     (y < prob_hi[1]) && (z > prob_lo[2]) && (z < prob_hi[2]));

                scale_fac_cc(i, j, k, 0) = in_domain ? fac_x : 1.0;
                scale_fac_cc(i, j, k, 1) = in_domain ? fac_y : 1.0;
                scale_fac_cc(i, j, k, 2) = in_domain ? fac_z : 1.0;

                scale_detJ_cc(i, j, k) = scale_fac_cc(i, j, k, 0) *
                                         scale_fac_cc(i, j, k, 1) *
                                         scale_fac_cc(i, j, k, 2);
            });

        const auto& nbx = mfi.grownnodaltilebox();
        amrex::Array4<amrex::Real> const& scale_fac_nd =
            (*m_mesh_scale_fac_nd)(lev).array(mfi);
        amrex::Array4<amrex::Real> const& scale_detJ_nd =
            (*m_mesh_scale_detJ_nd)(lev).array(mfi);
        amrex::ParallelFor(
            nbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                amrex::Real x = prob_lo[0] + i * dx[0];
                amrex::Real y = prob_lo[1] + j * dx[1];
                amrex::Real z = prob_lo[2] + k * dx[2];

                amrex::Real fac_x = eval_fac(x, beta[0], prob_lo[0], len[0]);
                amrex::Real fac_y = eval_fac(y, beta[1], prob_lo[1], len[1]);
                amrex::Real fac_z = eval_fac(z, beta[2], prob_lo[2], len[2]);

                bool in_domain =
                    ((x >= prob_lo[0] - eps) && (x <= prob_hi[0] + eps) &&
                     (y >= prob_lo[1] - eps) && (y <= prob_hi[1] + eps) &&
                     (z >= prob_lo[2] - eps) && (z <= prob_hi[2] + eps));

                scale_fac_nd(i, j, k, 0) = in_domain ? fac_x : 1.0;
                scale_fac_nd(i, j, k, 1) = in_domain ? fac_y : 1.0;
                scale_fac_nd(i, j, k, 2) = in_domain ? fac_z : 1.0;

                scale_detJ_nd(i, j, k) = scale_fac_nd(i, j, k, 0) *
                                         scale_fac_nd(i, j, k, 1) *
                                         scale_fac_nd(i, j, k, 2);
            });
    }

    // TODO: Call fill patch operators ?
}

/** Construct the mesh mapping field on cell faces
 */
void ChannelFlowMap::create_face_map(int lev, const amrex::Geometry& geom)
{
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> beta{
        {m_beta[0], m_beta[1], m_beta[2]}};
    const auto eps = m_eps;

    const auto& dx = geom.CellSizeArray();
    const auto& prob_lo = geom.ProbLoArray();
    const auto& prob_hi = geom.ProbHiArray();

    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> len{
        {prob_hi[0] - prob_lo[0], prob_hi[1] - prob_lo[1],
         prob_hi[2] - prob_lo[2]}};

    for (amrex::MFIter mfi((*m_mesh_scale_fac_xf)(lev)); mfi.isValid(); ++mfi) {

        const auto& bx = mfi.growntilebox();
        amrex::Array4<amrex::Real> const& scale_fac_xf =
            (*m_mesh_scale_fac_xf)(lev).array(mfi);
        amrex::Array4<amrex::Real> const& scale_detJ_xf =
            (*m_mesh_scale_detJ_xf)(lev).array(mfi);
        amrex::ParallelFor(
            bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                amrex::Real x = prob_lo[0] + i * dx[0];
                amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1];
                amrex::Real z = prob_lo[2] + (k + 0.5) * dx[2];

                amrex::Real fac_x = eval_fac(x, beta[0], prob_lo[0], len[0]);
                amrex::Real fac_y = eval_fac(y, beta[1], prob_lo[1], len[1]);
                amrex::Real fac_z = eval_fac(z, beta[2], prob_lo[2], len[2]);

                bool in_domain =
                    ((x >= prob_lo[0] - eps) && (x <= prob_hi[0] + eps) &&
                     (y > prob_lo[1]) && (y < prob_hi[1]) && (z > prob_lo[2]) &&
                     (z < prob_hi[2]));

                scale_fac_xf(i, j, k, 0) = in_domain ? fac_x : 1.0;
                scale_fac_xf(i, j, k, 1) = in_domain ? fac_y : 1.0;
                scale_fac_xf(i, j, k, 2) = in_domain ? fac_z : 1.0;

                scale_detJ_xf(i, j, k) = scale_fac_xf(i, j, k, 0) *
                                         scale_fac_xf(i, j, k, 1) *
                                         scale_fac_xf(i, j, k, 2);
            });
    }

    for (amrex::MFIter mfi((*m_mesh_scale_fac_yf)(lev)); mfi.isValid(); ++mfi) {

        const auto& bx = mfi.growntilebox();
        amrex::Array4<amrex::Real> const& scale_fac_yf =
            (*m_mesh_scale_fac_yf)(lev).array(mfi);
        amrex::Array4<amrex::Real> const& scale_detJ_yf =
            (*m_mesh_scale_detJ_yf)(lev).array(mfi);
        amrex::ParallelFor(
            bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];
                amrex::Real y = prob_lo[1] + j * dx[1];
                amrex::Real z = prob_lo[2] + (k + 0.5) * dx[2];

                amrex::Real fac_x = eval_fac(x, beta[0], prob_lo[0], len[0]);
                amrex::Real fac_y = eval_fac(y, beta[1], prob_lo[1], len[1]);
                amrex::Real fac_z = eval_fac(z, beta[2], prob_lo[2], len[2]);

                bool in_domain =
                    ((x > prob_lo[0]) && (x < prob_hi[0]) &&
                     (y >= prob_lo[1] - eps) && (y <= prob_hi[1] + eps) &&
                     (z > prob_lo[2]) && (z < prob_hi[2]));

                scale_fac_yf(i, j, k, 0) = in_domain ? fac_x : 1.0;
                scale_fac_yf(i, j, k, 1) = in_domain ? fac_y : 1.0;
                scale_fac_yf(i, j, k, 2) = in_domain ? fac_z : 1.0;

                scale_detJ_yf(i, j, k) = scale_fac_yf(i, j, k, 0) *
                                         scale_fac_yf(i, j, k, 1) *
                                         scale_fac_yf(i, j, k, 2);
            });
    }

    for (amrex::MFIter mfi((*m_mesh_scale_fac_zf)(lev)); mfi.isValid(); ++mfi) {

        const auto& bx = mfi.growntilebox();
        amrex::Array4<amrex::Real> const& scale_fac_zf =
            (*m_mesh_scale_fac_zf)(lev).array(mfi);
        amrex::Array4<amrex::Real> const& scale_detJ_zf =
            (*m_mesh_scale_detJ_zf)(lev).array(mfi);
        amrex::ParallelFor(
            bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];
                amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1];
                amrex::Real z = prob_lo[2] + k * dx[2];

                amrex::Real fac_x = eval_fac(x, beta[0], prob_lo[0], len[0]);
                amrex::Real fac_y = eval_fac(y, beta[1], prob_lo[1], len[1]);
                amrex::Real fac_z = eval_fac(z, beta[2], prob_lo[2], len[2]);

                bool in_domain =
                    ((x > prob_lo[0]) && (x < prob_hi[0]) && (y > prob_lo[1]) &&
                     (y < prob_hi[1]) && (z >= prob_lo[2] - eps) &&
                     (z <= prob_hi[2] + eps));

                scale_fac_zf(i, j, k, 0) = in_domain ? fac_x : 1.0;
                scale_fac_zf(i, j, k, 1) = in_domain ? fac_y : 1.0;
                scale_fac_zf(i, j, k, 2) = in_domain ? fac_z : 1.0;

                scale_detJ_zf(i, j, k) = scale_fac_zf(i, j, k, 0) *
                                         scale_fac_zf(i, j, k, 1) *
                                         scale_fac_zf(i, j, k, 2);
            });
    }

    // TODO: Call fill patch operators ?
}

/** Construct the non-uniform mesh field
 */
void ChannelFlowMap::create_non_uniform_mesh(
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

    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> beta{
        {m_beta[0], m_beta[1], m_beta[2]}};
    const auto eps = m_eps;

    const auto& dx = geom.CellSizeArray();
    const auto& prob_lo = geom.ProbLoArray();
    const auto& prob_hi = geom.ProbHiArray();

    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> len{
        {probhi_physical[0] - prob_lo[0], probhi_physical[1] - prob_lo[1],
         probhi_physical[2] - prob_lo[2]}};

    for (amrex::MFIter mfi((*m_non_uniform_coord_cc)(lev)); mfi.isValid();
         ++mfi) {

        const auto& bx = mfi.growntilebox();
        amrex::Array4<amrex::Real> const& nu_coord_cc =
            (*m_non_uniform_coord_cc)(lev).array(mfi);
        amrex::ParallelFor(
            bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];
                amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1];
                amrex::Real z = prob_lo[2] + (k + 0.5) * dx[2];

                amrex::Real x_non_uni =
                    eval_coord(x, beta[0], prob_lo[0], len[0]);
                amrex::Real y_non_uni =
                    eval_coord(y, beta[1], prob_lo[1], len[1]);
                amrex::Real z_non_uni =
                    eval_coord(z, beta[2], prob_lo[2], len[2]);

                bool in_domain =
                    ((x > prob_lo[0]) && (x < prob_hi[0]) && (y > prob_lo[1]) &&
                     (y < prob_hi[1]) && (z > prob_lo[2]) && (z < prob_hi[2]));

                nu_coord_cc(i, j, k, 0) = in_domain ? x_non_uni : x;
                nu_coord_cc(i, j, k, 1) = in_domain ? y_non_uni : y;
                nu_coord_cc(i, j, k, 2) = in_domain ? z_non_uni : z;
            });

        const auto& nbx = mfi.grownnodaltilebox();
        amrex::Array4<amrex::Real> const& nu_coord_nd =
            (*m_non_uniform_coord_nd)(lev).array(mfi);
        amrex::ParallelFor(
            nbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                amrex::Real x = prob_lo[0] + i * dx[0];
                amrex::Real y = prob_lo[1] + j * dx[1];
                amrex::Real z = prob_lo[2] + k * dx[2];

                amrex::Real x_non_uni =
                    eval_coord(x, beta[0], prob_lo[0], len[0]);
                amrex::Real y_non_uni =
                    eval_coord(y, beta[1], prob_lo[1], len[1]);
                amrex::Real z_non_uni =
                    eval_coord(z, beta[2], prob_lo[2], len[2]);

                bool in_domain =
                    ((x >= prob_lo[0] - eps) && (x <= prob_hi[0] + eps) &&
                     (y >= prob_lo[1] - eps) && (y <= prob_hi[1] + eps) &&
                     (z >= prob_lo[2] - eps) && (z <= prob_hi[2] + eps));

                nu_coord_nd(i, j, k, 0) = in_domain ? x_non_uni : x;
                nu_coord_nd(i, j, k, 1) = in_domain ? y_non_uni : y;
                nu_coord_nd(i, j, k, 2) = in_domain ? z_non_uni : z;
            });
    }
}

} // namespace amr_wind::channel_map
