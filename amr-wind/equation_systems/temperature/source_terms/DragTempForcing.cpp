#include "amr-wind/equation_systems/temperature/source_terms/DragTempForcing.H"
#include "amr-wind/utilities/IOManager.H"
#include "amr-wind/wind_energy/MOData.H"
#include "AMReX_ParmParse.H"
#include "AMReX_Gpu.H"
#include "AMReX_Random.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

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
    pp.query("bc_forcing_time_factor", m_forcing_time_factor);
    amrex::ParmParse pp_abl("ABL");
    pp_abl.query("wall_het_model", m_wall_het_model);
    pp_abl.query("monin_obukhov_length", m_monin_obukhov_length);
    pp_abl.query("kappa", m_kappa);
    pp_abl.query("mo_gamma_m", m_gamma_m);
    pp_abl.query("mo_beta_m", m_beta_m);
    pp_abl.query("mo_gamma_m", m_gamma_h);
    pp_abl.query("mo_beta_m", m_beta_h);

    amrex::ParmParse pp_incflo("incflo");
    pp_incflo.queryarr("gravity", m_gravity);
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
    const amrex::Real z0_min = 1.0e-4_rt;
    const amrex::Real monin_obukhov_length = m_monin_obukhov_length;
    const auto& dt = m_time.delta_t();
    const bool is_mol = (m_wall_het_model == "mol");
    const amrex::Real beta_m = m_beta_m;
    const amrex::Real gamma_m = m_gamma_m;
    const amrex::Real beta_h = m_beta_h;
    const amrex::Real gamma_h = m_gamma_h;
    // Displacement height (e.g., tree canopy); log law uses z - d
    const bool has_displacement =
        this->m_sim.repo().field_exists("terrain_displacement");
    const auto& displacement = has_displacement
                                   ? this->m_sim.repo()
                                         .get_field("terrain_displacement")(lev)
                                         .const_array(mfi)
                                   : amrex::Array4<amrex::Real const>();
    const auto tiny = std::numeric_limits<amrex::Real>::epsilon();
    const amrex::Real cd_max = 10.0_rt;
    const amrex::Real T0 = m_soil_temperature;
    const amrex::Real time_factor = m_forcing_time_factor;
    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        const amrex::Real z0 =
            amrex::max<amrex::Real>(terrainz0(i, j, k), z0_min);
        const amrex::Real disp =
            has_displacement ? displacement(i, j, k) : 0.0_rt;
        const amrex::Real z_neighbour =
            MOData::effective_height(1.5_rt * dx[2], disp, z0);
        const amrex::Real z_cell =
            MOData::effective_height(0.5_rt * dx[2], disp, z0);
        const amrex::Real psi_m =
            is_mol ? MOData::calc_psi_m(
                         z_neighbour / monin_obukhov_length, beta_m, gamma_m)
                   : 0.0_rt;
        const amrex::Real psi_h_neighbour =
            is_mol ? MOData::calc_psi_h(
                         z_neighbour / monin_obukhov_length, beta_h, gamma_h)
                   : 0.0_rt;
        const amrex::Real psi_h_cell =
            is_mol ? MOData::calc_psi_h(
                         z_cell / monin_obukhov_length, beta_h, gamma_h)
                   : 0.0_rt;
        const amrex::Real ux1 = vel(i, j, k, 0);
        const amrex::Real uy1 = vel(i, j, k, 1);
        const amrex::Real uz1 = vel(i, j, k, 2);
        const amrex::Real theta = temperature(i, j, k, 0);
        const amrex::Real theta2 = temperature(i, j, k + 1, 0);
        const amrex::Real wspd = std::sqrt((ux1 * ux1) + (uy1 * uy1));
        const amrex::Real ustar =
            wspd * kappa / (std::log(z_neighbour / z0) - psi_m);
        //! We do not know the actual temperature so use cell above
        const amrex::Real thetastar =
            theta * ustar * ustar /
            (kappa * gravity_mod * monin_obukhov_length);
        const amrex::Real surf_temp =
            theta2 - (thetastar / kappa *
                      (std::log(z_neighbour / z0) - psi_h_neighbour));
        const amrex::Real tTarget =
            surf_temp +
            (thetastar / kappa * (std::log(z_cell / z0) - psi_h_cell));
        const amrex::Real bc_forcing_t =
            -(tTarget - theta) / (time_factor * dt);
        const amrex::Real m =
            std::sqrt((ux1 * ux1) + (uy1 * uy1) + (uz1 * uz1));
        const amrex::Real Cd = amrex::min<amrex::Real>(
            drag_coefficient / (m + tiny), cd_max / dx[2]);
        src_term(i, j, k, 0) -=
            ((Cd * (theta - T0) * blank(i, j, k, 0)) +
             (bc_forcing_t * drag(i, j, k)));
    });
}

} // namespace amr_wind::pde::temperature
