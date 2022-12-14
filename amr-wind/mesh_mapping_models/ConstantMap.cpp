#include "amr-wind/mesh_mapping_models/ConstantMap.H"

#include "AMReX_ParmParse.H"

namespace amr_wind::const_map {

ConstantMap::ConstantMap()
{
    amrex::ParmParse pp("ConstantMap");
    pp.queryarr("scaling_factor", m_fac, 0, AMREX_SPACEDIM);
}

/** Construct the mesh mapping field
 */
void ConstantMap::create_map(int lev, const amrex::Geometry& geom)
{
    create_cell_node_map(lev);
    create_face_map(lev);
    create_non_uniform_mesh(lev, geom);
}

/** Construct the mesh mapping field on cell centers and nodes
 */
void ConstantMap::create_cell_node_map(int lev)
{
    amrex::Real fac_x = m_fac[0];
    amrex::Real fac_y = m_fac[1];
    amrex::Real fac_z = m_fac[2];

    for (amrex::MFIter mfi((*m_mesh_scale_fac_cc)(lev)); mfi.isValid(); ++mfi) {

        const auto& bx = mfi.growntilebox();
        amrex::Array4<amrex::Real> const& scale_fac_cc =
            (*m_mesh_scale_fac_cc)(lev).array(mfi);
        amrex::Array4<amrex::Real> const& scale_detJ_cc =
            (*m_mesh_scale_detJ_cc)(lev).array(mfi);
        amrex::ParallelFor(
            bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                scale_fac_cc(i, j, k, 0) = fac_x;
                scale_fac_cc(i, j, k, 1) = fac_y;
                scale_fac_cc(i, j, k, 2) = fac_z;
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
                scale_fac_nd(i, j, k, 0) = fac_x;
                scale_fac_nd(i, j, k, 1) = fac_y;
                scale_fac_nd(i, j, k, 2) = fac_z;
                scale_detJ_nd(i, j, k) = scale_fac_nd(i, j, k, 0) *
                                         scale_fac_nd(i, j, k, 1) *
                                         scale_fac_nd(i, j, k, 2);
            });
    }
}

/** Construct the mesh mapping field on cell faces
 */
void ConstantMap::create_face_map(int lev)
{
    amrex::Real fac_x = m_fac[0];
    amrex::Real fac_y = m_fac[1];
    amrex::Real fac_z = m_fac[2];

    for (amrex::MFIter mfi((*m_mesh_scale_fac_xf)(lev)); mfi.isValid(); ++mfi) {
        const auto& bx = mfi.growntilebox();
        amrex::Array4<amrex::Real> const& scale_fac_xf =
            (*m_mesh_scale_fac_xf)(lev).array(mfi);
        amrex::Array4<amrex::Real> const& scale_detJ_xf =
            (*m_mesh_scale_detJ_xf)(lev).array(mfi);
        amrex::ParallelFor(
            bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                scale_fac_xf(i, j, k, 0) = fac_x;
                scale_fac_xf(i, j, k, 1) = fac_y;
                scale_fac_xf(i, j, k, 2) = fac_z;
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
                scale_fac_yf(i, j, k, 0) = fac_x;
                scale_fac_yf(i, j, k, 1) = fac_y;
                scale_fac_yf(i, j, k, 2) = fac_z;
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
                scale_fac_zf(i, j, k, 0) = fac_x;
                scale_fac_zf(i, j, k, 1) = fac_y;
                scale_fac_zf(i, j, k, 2) = fac_z;
                scale_detJ_zf(i, j, k) = scale_fac_zf(i, j, k, 0) *
                                         scale_fac_zf(i, j, k, 1) *
                                         scale_fac_zf(i, j, k, 2);
            });
    }
}

/** Construct the non-uniform mesh field
 */
void ConstantMap::create_non_uniform_mesh(int lev, const amrex::Geometry& geom)
{

    const auto& problo = geom.ProbLoArray();
    const auto& dx = geom.CellSizeArray();

    for (amrex::MFIter mfi((*m_non_uniform_coord_cc)(lev)); mfi.isValid();
         ++mfi) {

        const auto& bx = mfi.growntilebox();
        amrex::Array4<amrex::Real> const& scale_fac_cc =
            (*m_mesh_scale_fac_cc)(lev).array(mfi);
        amrex::Array4<amrex::Real> const& nu_coord_cc =
            (*m_non_uniform_coord_cc)(lev).array(mfi);
        amrex::ParallelFor(
            bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                nu_coord_cc(i, j, k, 0) =
                    problo[0] + (i + 0.5) * dx[0] * scale_fac_cc(i, j, k, 0);
                nu_coord_cc(i, j, k, 1) =
                    problo[1] + (j + 0.5) * dx[1] * scale_fac_cc(i, j, k, 1);
                nu_coord_cc(i, j, k, 2) =
                    problo[2] + (k + 0.5) * dx[2] * scale_fac_cc(i, j, k, 2);
            });

        const auto& nbx = mfi.grownnodaltilebox();
        amrex::Array4<amrex::Real> const& scale_fac_nd =
            (*m_mesh_scale_fac_nd)(lev).array(mfi);
        amrex::Array4<amrex::Real> const& nu_coord_nd =
            (*m_non_uniform_coord_nd)(lev).array(mfi);
        amrex::ParallelFor(
            nbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                nu_coord_nd(i, j, k, 0) =
                    problo[0] + i * dx[0] * scale_fac_nd(i, j, k, 0);
                nu_coord_nd(i, j, k, 1) =
                    problo[1] + j * dx[1] * scale_fac_nd(i, j, k, 1);
                nu_coord_nd(i, j, k, 2) =
                    problo[2] + k * dx[2] * scale_fac_nd(i, j, k, 2);
            });
    }
}

} // namespace amr_wind::const_map
