
#include "amr-wind/equation_systems/temperature/source_terms/DragTempForcing.H"
#include "amr-wind/utilities/IOManager.H"

#include "AMReX_ParmParse.H"
#include "AMReX_Gpu.H"
#include "AMReX_Random.H"

namespace amr_wind::pde::temperature {

DragTempForcing::DragTempForcing(const CFDSim& sim)
    : m_time(sim.time())
    , m_sim(sim)
    , m_mesh(sim.mesh())
    , m_velocity(sim.repo().get_field("velocity"))
    , m_temperature(sim.repo().get_field("temperature"))
{
    amrex::ParmParse pp("DragTempForcing");
    pp.query("drag_coefficient", m_drag_coefficient);
    pp.query("internal_temperature", m_internal_temperature);
    amrex::ParmParse pp_abl("ABL");
    pp_abl.query("wall_het_model", m_wall_het_model);
    pp_abl.query("mol_length", m_mol_length);
    pp_abl.query("surface_roughness_z0", m_z0);
    pp_abl.query("kappa", m_kappa);
    pp_abl.query("mo_gamma_m", m_gamma_m);
    pp_abl.query("mo_beta_m", m_beta_m);
    pp_abl.query("mo_gamma_m", m_gamma_h);
    pp_abl.query("mo_beta_m", m_beta_h);
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
    auto* const m_terrain_drag =
        &this->m_sim.repo().get_int_field("terrain_drag");
    const auto& drag = (*m_terrain_drag)(lev).const_array(mfi);
    const auto& geom = m_mesh.Geom(lev);
    const auto& dx = geom.CellSizeArray();
    const amrex::Real drag_coefficient = m_drag_coefficient / dx[2];
    const amrex::Real internal_temperature = m_internal_temperature;
    const amrex::Real gravity_mod = 9.81;
    amrex::Real psi_m = 0.0;
    amrex::Real psi_h_neighbour = 0.0;
    amrex::Real psi_h_cell = 0.0;
    const amrex::Real kappa = m_kappa;
    const amrex::Real z0 = m_z0;
    const amrex::Real mol_length = m_mol_length;
    const auto& dt = m_time.delta_t();
    if (m_wall_het_model == "mol") {
        psi_m = stability(1.5 * dx[2], m_mol_length, m_gamma_m, m_beta_m);
        psi_h_neighbour =
            thermal_stability(1.5 * dx[2], m_mol_length, m_gamma_h, m_beta_h);
        psi_h_cell =
            thermal_stability(0.5 * dx[2], m_mol_length, m_gamma_h, m_beta_h);
    }
    const auto tiny = std::numeric_limits<amrex::Real>::epsilon();
    const amrex::Real cd_max = 10.0;
    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        const amrex::Real ux1 = vel(i, j, k, 0);
        const amrex::Real uy1 = vel(i, j, k, 1);
        const amrex::Real uz1 = vel(i, j, k, 2);
        const amrex::Real theta = temperature(i, j, k, 0);
        const amrex::Real theta2 = temperature(i, j, k + 1, 0);
        const amrex::Real wspd = std::sqrt(ux1 * ux1 + uy1 * uy1);
        const amrex::Real ustar =
            wspd * kappa / (std::log(1.5 * dx[2] / z0) - psi_m);
        //! We do not know the actual temperature so use cell above
        const amrex::Real thetastar =
            theta * ustar * ustar / (kappa * gravity_mod * mol_length);
        const amrex::Real surf_temp =
            theta2 -
            thetastar / kappa * (std::log(1.5 * dx[2] / z0) - psi_h_neighbour);
        const amrex::Real tTarget =
            surf_temp +
            thetastar / kappa * (std::log(0.5 * dx[2] / z0) - psi_h_cell);
        const amrex::Real bc_forcing_t = -(tTarget - theta) / dt;
        const amrex::Real m = std::sqrt(ux1 * ux1 + uy1 * uy1 + uz1 * uz1);
        const amrex::Real Cd =
            std::min(drag_coefficient / (m + tiny), cd_max / dx[2]);
        src_term(i, j, k, 0) -=
            (Cd * (theta - internal_temperature) * blank(i, j, k, 0) +
             bc_forcing_t * drag(i, j, k));
    });
}

} // namespace amr_wind::pde::temperature
