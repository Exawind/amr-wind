
#include "amr-wind/equation_systems/temperature/source_terms/DragTempForcing.H"
#include "amr-wind/utilities/IOManager.H"

#include "AMReX_ParmParse.H"
#include "AMReX_Gpu.H"
#include "AMReX_Random.H"

namespace amr_wind::pde::temperature {

DragTempForcing::DragTempForcing(const CFDSim& sim)
    : m_sim(sim)
    , m_mesh(sim.mesh())
    , m_velocity(sim.repo().get_field("velocity"))
    , m_temperature(sim.repo().get_field("temperature"))
{
    amrex::ParmParse pp("DragTempForcing");
    pp.query("dragCoefficient", m_drag);
    pp.query("RefT", m_internalRefT);
}

DragTempForcing::~DragTempForcing() = default;

void DragTempForcing::operator()(
    const int lev,
    const amrex::MFIter& mfi,
    const amrex::Box& bx,
    const FieldState fstate,
    const amrex::Array4<amrex::Real>& src_term) const
{
    const auto& vel =
        m_velocity.state(field_impl::dof_state(fstate))(lev).const_array(mfi);
    const auto& temperature =
        m_temperature.state(field_impl::dof_state(fstate))(lev).const_array(
            mfi);
    const bool is_terrain =
        this->m_sim.repo().int_field_exists("terrain_blank");
    if (!is_terrain) {
        amrex::Abort("Need terrain blanking variable to use this source term");
    }
    auto* const m_terrain_blank =
        &this->m_sim.repo().get_int_field("terrain_blank");
    const auto& blank = (*m_terrain_blank)(lev).const_array(mfi);
    const auto& geom = m_mesh.Geom(lev);
    const auto& dx = geom.CellSizeArray();
    const amrex::Real gpu_drag = m_drag / dx[2];
    const amrex::Real gpu_TRef = m_internalRefT;
    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        const amrex::Real ux1 = vel(i, j, k, 0);
        const amrex::Real uy1 = vel(i, j, k, 1);
        const amrex::Real uz1 = vel(i, j, k, 2);
        const amrex::Real m = std::sqrt(ux1 * ux1 + uy1 * uy1 + uz1 * uz1);
        const amrex::Real Cd = std::min(gpu_drag / (m + 1e-5), 10 / dx[2]);
        src_term(i, j, k, 0) -=
            (Cd * (temperature(i, j, k, 0) - gpu_TRef) * blank(i, j, k, 0));
    });
}

} // namespace amr_wind::pde::temperature
