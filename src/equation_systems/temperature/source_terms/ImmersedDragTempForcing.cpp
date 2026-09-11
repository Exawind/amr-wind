#include "src/equation_systems/temperature/source_terms/ImmersedDragTempForcing.H"
#include "src/equation_systems/icns/source_terms/ImmersedDragForcing.H"
#include "src/physics/ImmersedTerrain.H"
#include "src/physics/ImmersedWallModel.H"
#include "src/utilities/constants.H"
#include "AMReX_Gpu.H"
#include "AMReX_ParmParse.H"
#include "AMReX_REAL.H"

#include <cstdint>

using namespace amrex::literals;
using kynema_sgf::immersed_wall::WallModel;
using kynema_sgf::immersed_wall::WallParams;
using kynema_sgf::immersed_wall::WallPatch;
using kynema_sgf::immersedterrain::ImmersedTerrain;

namespace {
enum class SurfaceCondition : std::uint8_t {
    obukhov_length = 0,
    surface_temperature,
    heat_flux
};
} // namespace

namespace kynema_sgf::pde::temperature {

ImmersedDragTempForcing::ImmersedDragTempForcing(const CFDSim& sim)
    : m_time(sim.time())
    , m_sim(sim)
    , m_mesh(sim.mesh())
    , m_velocity(sim.repo().get_field("velocity"))
    , m_temperature(sim.repo().get_field("temperature"))
{
    amrex::ParmParse pp(identifier());
    pp.query("drag_coefficient", m_drag_coefficient);
    pp.query("soil_temperature", m_soil_temperature);
    pp.query("surface_condition", m_surface_condition);
    pp.query("surface_heat_flux", m_surface_heat_flux);
    pp.query("force_laminar", m_force_laminar);
    if (m_surface_condition != "obukhov_length" &&
        m_surface_condition != "surface_temperature" &&
        m_surface_condition != "heat_flux") {
        amrex::Abort(
            identifier() +
            ".surface_condition must be obukhov_length, "
            "surface_temperature or heat_flux");
    }

    // Wall-model geometry and time scale are owned by the momentum source
    amrex::ParmParse pp_mom(icns::ImmersedDragForcing::identifier());
    pp_mom.query("wall_model", m_wall_model);
    pp_mom.query("reference_distance", m_reference_distance);
    pp_mom.query("bc_forcing_time_scale", m_bc_forcing_time_scale);
    pp_mom.query("bc_forcing_time_factor", m_forcing_time_factor);
    pp_mom.query("minimum_z0", m_min_z0);
    amrex::ParmParse pp_terrain(ImmersedTerrain::identifier());
    pp_terrain.query("solid_threshold", m_solid_threshold);
    pp_terrain.query("drag_weight", m_drag_weight);

    std::string turbulence_model = "Laminar";
    amrex::ParmParse pp_turb("turbulence");
    pp_turb.query("model", turbulence_model);
    m_is_laminar = (turbulence_model == "Laminar") || m_force_laminar;

    amrex::ParmParse pp_abl("ABL");
    pp_abl.query("wall_het_model", m_wall_het_model);
    pp_abl.query("monin_obukhov_length", m_monin_obukhov_length);
    pp_abl.query("kappa", m_kappa);
    pp_abl.query("mo_gamma_m", m_gamma_m);
    pp_abl.query("mo_beta_m", m_beta_m);
    pp_abl.query("mo_gamma_h", m_gamma_h);
    pp_abl.query("mo_beta_h", m_beta_h);

    amrex::ParmParse pp_incflo("incflo");
    pp_incflo.queryarr("gravity", m_gravity);

    if (!sim.repo().field_exists("terrain_fraction")) {
        amrex::Abort(
            identifier() +
            " requires the ImmersedTerrain physics "
            "(terrain_fraction field not found)");
    }
}

ImmersedDragTempForcing::~ImmersedDragTempForcing() = default;

void ImmersedDragTempForcing::operator()(
    const int lev, const FieldState fstate, amrex::MultiFab& src_term) const
{
    const auto& repo = m_sim.repo();
    auto const& src_arrs = src_term.arrays();
    auto const& vel_arrs =
        m_velocity.state(field_impl::dof_state(fstate))(lev).const_arrays();
    auto const& temp_arrs =
        m_temperature.state(field_impl::dof_state(fstate))(lev).const_arrays();
    auto const& frac_arrs =
        repo.get_field("terrain_fraction")(lev).const_arrays();
    auto const& mask_arrs =
        repo.get_int_field("terrain_mask")(lev).const_arrays();
    auto const& surf_arrs =
        repo.get_field("terrain_surface")(lev).const_arrays();
    auto const& z0_arrs =
        repo.get_field("terrain_roughness")(lev).const_arrays();

    const auto& geom = m_mesh.Geom(lev);
    const auto dx = geom.CellSizeArray();
    const auto prob_lo = geom.ProbLoArray();
    const amrex::Real dt = m_time.delta_t();

    const amrex::Real relax_rate = m_drag_coefficient / dx[2];
    const amrex::Real theta_soil = m_soil_temperature;
    const amrex::Real min_z0 = m_min_z0;
    const amrex::Real solid_threshold = m_solid_threshold;
    const bool apply_wall_model = !m_is_laminar;
    const WallModel wall_model =
        kynema_sgf::immersed_wall::parse_wall_model(m_wall_model);
    const int actual_reference = (m_reference_distance == "actual") ? 1 : 0;
    const int center_weight = (m_drag_weight == "center") ? 1 : 0;
    const SurfaceCondition condition =
        (m_surface_condition == "surface_temperature")
            ? SurfaceCondition::surface_temperature
        : (m_surface_condition == "heat_flux")
            ? SurfaceCondition::heat_flux
            : SurfaceCondition::obukhov_length;
    const amrex::Real surface_heat_flux = m_surface_heat_flux;
    const amrex::Real gravity = std::abs(m_gravity[2]);
    // theta* from a prescribed L: theta u*^2 / (kappa g L)
    const amrex::Real inv_kappa_g_L =
        1.0_rt / (m_kappa * gravity * m_monin_obukhov_length);

    WallParams wp{};
    wp.kappa = m_kappa;
    wp.beta_m = m_beta_m;
    wp.gamma_m = m_gamma_m;
    wp.beta_h = m_beta_h;
    wp.gamma_h = m_gamma_h;
    wp.inv_mo_length =
        (m_wall_het_model == "mol") ? 1.0_rt / m_monin_obukhov_length : 0.0_rt;
    wp.dt = dt;
    wp.tau_floor = m_forcing_time_factor * dt;
    wp.use_wall_time_scale = (m_bc_forcing_time_scale == "wall");

    amrex::ParallelFor(
        src_term, [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
            const auto& vel = vel_arrs[nbx];
            const auto& temp = temp_arrs[nbx];
            const auto& frac = frac_arrs[nbx];
            const auto& src = src_arrs[nbx];

            const amrex::Real beta = frac(i, j, k, 0);
            const amrex::Real w_solid = kynema_sgf::immersed_wall::solid_weight(
                beta, center_weight != 0, solid_threshold);
            const int cell_mask = mask_arrs[nbx](i, j, k, 0);
            const amrex::Real theta = temp(i, j, k, 0);

            // 1. Relaxation toward the soil temperature inside the terrain
            if (w_solid > 0.0_rt) {
                const amrex::Real C_eff =
                    kynema_sgf::immersed_wall::exact_relaxation_rate(
                        w_solid * relax_rate, dt);
                src(i, j, k, 0) -= C_eff * (theta - theta_soil);
            }

            // 2. Heat-flux wall model in surface cells
            if (!apply_wall_model ||
                cell_mask != ImmersedTerrain::mask_surface ||
                w_solid >= 1.0_rt) {
                return;
            }
            const amrex::Real z0 =
                amrex::max<amrex::Real>(z0_arrs[nbx](i, j, k, 0), min_z0);
            const amrex::Real z_c = prob_lo[2] + ((k + 0.5_rt) * dx[2]);

            amrex::GpuArray<WallPatch, 2 * AMREX_SPACEDIM> patches{};
            const int np = kynema_sgf::immersed_wall::wall_patches(
                wall_model, i, j, k, beta, frac, surf_arrs[nbx], dx, z_c, z0,
                solid_threshold, actual_reference != 0, patches.data());

            amrex::Real force = 0.0_rt;
            amrex::Real weight_sum = 0.0_rt;
            for (int ip = 0; ip < np; ++ip) {
                const WallPatch& p = patches[ip];
                const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> u_ref{
                    vel(p.ir, p.jr, p.kr, 0), vel(p.ir, p.jr, p.kr, 1),
                    vel(p.ir, p.jr, p.kr, 2)};
                const auto ut_ref =
                    kynema_sgf::immersed_wall::tangential(u_ref, p.nrm);
                const amrex::Real ustar =
                    kynema_sgf::immersed_wall::friction_velocity(
                        kynema_sgf::immersed_wall::magnitude(ut_ref), p, z0,
                        wp);
                const amrex::Real theta_ref = temp(p.ir, p.jr, p.kr, 0);
                const amrex::Real dd1 = amrex::max<amrex::Real>(p.d1, z0);
                const amrex::Real dd2 = amrex::max<amrex::Real>(p.d2, dd1);
                const amrex::Real phi_h2 = wp.phi_h(dd2, z0);

                amrex::Real theta_star = 0.0_rt;
                amrex::Real theta_s = theta_soil;
                if (condition == SurfaceCondition::surface_temperature) {
                    theta_star = wp.kappa * (theta_ref - theta_soil) / phi_h2;
                } else if (condition == SurfaceCondition::heat_flux) {
                    theta_star = -surface_heat_flux /
                                 (ustar + kynema_sgf::constants::EPS);
                    theta_s = theta_ref - (theta_star / wp.kappa * phi_h2);
                } else {
                    theta_star = theta * ustar * ustar * inv_kappa_g_L;
                    theta_s = theta_ref - (theta_star / wp.kappa * phi_h2);
                }
                force += p.weight *
                         kynema_sgf::immersed_wall::wall_heat_forcing(
                             theta, theta_s, theta_star, ustar, p, z0, wp);
                weight_sum += p.weight;
            }
            if (weight_sum > 0.0_rt) {
                src(i, j, k, 0) += (1.0_rt - w_solid) * force / weight_sum;
            }
        });
    amrex::Gpu::streamSynchronize();
}

} // namespace kynema_sgf::pde::temperature
