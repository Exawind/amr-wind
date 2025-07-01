
#include "amr-wind/equation_systems/temperature/source_terms/DragTempForcing.H"
#include "amr-wind/utilities/IOManager.H"
#include "amr-wind/wind_energy/MOData.H"

#include "AMReX_ParmParse.H"
#include "AMReX_Gpu.H"
#include "AMReX_Random.H"
#include "amr-wind/utilities/constants.H"

namespace {
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real compute_target_theta(
    const amrex::Real ux,
    const amrex::Real uy,
    const amrex::Real theta,
    const amrex::Real monin_obukhov_length,
    const amrex::Real gravity_mod,
    const amrex::Real dx,
    const amrex::Real z0,
    const amrex::Real kappa)
{
    const amrex::Real wspd = std::sqrt(ux * ux + uy * uy);
    const amrex::Real ustar = wspd * kappa / std::log(1.5 * dx / z0);
    const amrex::Real thetastar =
        theta * ustar * ustar / (kappa * gravity_mod * monin_obukhov_length);
    const amrex::Real surf_temp =
        theta - thetastar / kappa * (std::log(1.5 * dx / z0));
    const amrex::Real tTarget =
        surf_temp + thetastar / kappa * (std::log(0.5 * dx / z0));
    return tTarget;
}
} // namespace
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
    pp.query("soil_temperature", m_soil_temperature);
    amrex::ParmParse pp_abl("ABL");
    pp_abl.query("wall_het_model", m_wall_het_model);
    pp_abl.query("monin_obukhov_length", m_monin_obukhov_length);
    pp_abl.query("kappa", m_kappa);
    pp_abl.query("mo_gamma_m", m_gamma_m);
    pp_abl.query("mo_beta_m", m_beta_m);
    pp_abl.query("mo_gamma_m", m_gamma_h);
    pp_abl.query("mo_beta_m", m_beta_h);
    {
        amrex::ParmParse pp_incflow("incflo");
        pp_incflow.queryarr("gravity", m_gravity);
    }
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
    auto* const m_terrainz0 = &this->m_sim.repo().get_field("terrainz0");
    const auto& terrainz0 = (*m_terrainz0)(lev).const_array(mfi);
    const auto& geom = m_mesh.Geom(lev);
    const auto& dx = geom.CellSizeArray();
    const amrex::Real drag_coefficient = m_drag_coefficient / dx[2];
    const amrex::Real gravity_mod = std::abs(m_gravity[2]);
    const amrex::Real kappa = m_kappa;
    const amrex::Real z0_min = 1e-4;
    const amrex::Real monin_obukhov_length = m_monin_obukhov_length;
    const auto& dt = m_time.delta_t();
    const amrex::Real psi_m =
        (m_wall_het_model == "mol")
            ? MOData::calc_psi_m(
                  1.5 * dx[2] / m_monin_obukhov_length, m_beta_m, m_gamma_m)
            : 0.0;
    const amrex::Real psi_h_neighbour =
        (m_wall_het_model == "mol")
            ? MOData::calc_psi_h(
                  1.5 * dx[2] / m_monin_obukhov_length, m_beta_h, m_gamma_h)
            : 0.0;
    const amrex::Real psi_h_cell =
        (m_wall_het_model == "mol")
            ? MOData::calc_psi_h(
                  0.5 * dx[2] / m_monin_obukhov_length, m_beta_h, m_gamma_h)
            : 0.0;
    const auto tiny = std::numeric_limits<amrex::Real>::epsilon();
    const amrex::Real cd_max = 10.0;
    const amrex::Real T0 = m_soil_temperature;
    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        const int cell_drag = (drag(i, j, k) > 0) ? 1 : 0;
        const amrex::Real z0 = std::max(terrainz0(i, j, k), z0_min);
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
            theta * ustar * ustar /
            (kappa * gravity_mod * monin_obukhov_length);
        const amrex::Real surf_temp =
            theta2 -
            thetastar / kappa * (std::log(1.5 * dx[2] / z0) - psi_h_neighbour);
        const amrex::Real tTarget =
            surf_temp +
            thetastar / kappa * (std::log(0.5 * dx[2] / z0) - psi_h_cell);
        amrex::Real bc_forcing_t = -(tTarget - theta) / dt;
        if (drag(i, j, k) > 1) {
            //! West
            amrex::Real tmp_temp_target = compute_target_theta(
                vel(i - 1, j, k, 2), vel(i - 1, j, k, 1), theta,
                monin_obukhov_length, gravity_mod, dx[0], z0, kappa);
            bc_forcing_t +=
                -(tmp_temp_target - theta) / dt * blank(i - 1, j, k);
            //! East
            tmp_temp_target = compute_target_theta(
                vel(i + 1, j, k, 2), vel(i + 1, j, k, 1), theta,
                monin_obukhov_length, gravity_mod, dx[0], z0, kappa);
            bc_forcing_t +=
                -(tmp_temp_target - theta) / dt * blank(i + 1, j, k);
            //! South
            tmp_temp_target = compute_target_theta(
                vel(i, j - 1, k, 2), vel(i, j - 1, k, 0), theta,
                monin_obukhov_length, gravity_mod, dx[1], z0, kappa);
            bc_forcing_t +=
                -(tmp_temp_target - theta) / dt * blank(i, j - 1, k);
            //! North
            tmp_temp_target = compute_target_theta(
                vel(i, j + 1, k, 2), vel(i, j + 1, k, 0), theta,
                monin_obukhov_length, gravity_mod, dx[1], z0, kappa);
            bc_forcing_t +=
                -(tmp_temp_target - theta) / dt * blank(i, j + 1, k);
        }
        const amrex::Real m = std::sqrt(ux1 * ux1 + uy1 * uy1 + uz1 * uz1);
        const amrex::Real Cd =
            std::min(drag_coefficient / (m + tiny), cd_max / dx[2]);
        src_term(i, j, k, 0) -=
            (Cd * (theta - T0) * blank(i, j, k, 0) + bc_forcing_t * cell_drag);
    });
}

} // namespace amr_wind::pde::temperature
