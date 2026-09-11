#include "src/equation_systems/icns/source_terms/ImmersedDragForcing.H"
#include "src/physics/ImmersedTerrain.H"
#include "src/physics/ImmersedWallModel.H"
#include "src/utilities/constants.H"
#include "AMReX_Gpu.H"
#include "AMReX_ParmParse.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;
using kynema_sgf::immersed_wall::WallModel;
using kynema_sgf::immersed_wall::WallParams;
using kynema_sgf::immersed_wall::WallPatch;
using kynema_sgf::immersedterrain::ImmersedTerrain;

namespace kynema_sgf::pde::icns {

ImmersedDragForcing::ImmersedDragForcing(const CFDSim& sim)
    : m_time(sim.time())
    , m_sim(sim)
    , m_mesh(sim.mesh())
    , m_velocity(sim.repo().get_field("velocity"))
{
    amrex::ParmParse pp(identifier());
    pp.query("drag_coefficient", m_drag_coefficient);
    pp.query("bc_forcing_time_factor", m_forcing_time_factor);
    pp.query("bc_forcing_time_scale", m_bc_forcing_time_scale);
    pp.query("wall_model", m_wall_model);
    pp.query("reference_distance", m_reference_distance);
    pp.query("minimum_z0", m_min_z0);
    if (m_reference_distance != "nominal" && m_reference_distance != "actual") {
        amrex::Abort(
            identifier() + ".reference_distance must be nominal or actual");
    }
    pp.query("force_laminar", m_force_laminar);
    if (m_bc_forcing_time_scale != "wall" &&
        m_bc_forcing_time_scale != "time_step") {
        amrex::Abort(
            identifier() + ".bc_forcing_time_scale must be wall or time_step");
    }
    if (m_wall_model != "cell_offset" && m_wall_model != "terrain_height" &&
        m_wall_model != "surface_normal") {
        amrex::Abort(
            identifier() +
            ".wall_model must be cell_offset, terrain_height or "
            "surface_normal");
    }

    // Same threshold the terrain physics uses to classify surface cells
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

    if (!sim.repo().field_exists("terrain_fraction")) {
        amrex::Abort(
            identifier() +
            " requires the ImmersedTerrain physics "
            "(terrain_fraction field not found)");
    }
}

ImmersedDragForcing::~ImmersedDragForcing() = default;

void ImmersedDragForcing::operator()(
    const int lev, const FieldState fstate, amrex::MultiFab& src_term) const
{
    const auto& repo = m_sim.repo();
    auto const& src_arrs = src_term.arrays();
    auto const& vel_arrs =
        m_velocity.state(field_impl::dof_state(fstate))(lev).const_arrays();
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

    // With ImmersedTerrain.implicit_projection the drag is applied inside the
    // projections through terrain_drag_rate; only the wall model remains here
    const bool implicit_drag = repo.field_exists("terrain_drag_rate");
    const amrex::Real drag_rate =
        implicit_drag ? 0.0_rt : m_drag_coefficient / dx[2];
    const amrex::Real min_z0 = m_min_z0;
    const amrex::Real solid_threshold = m_solid_threshold;
    const bool apply_wall_model = !m_is_laminar;
    const WallModel wall_model =
        kynema_sgf::immersed_wall::parse_wall_model(m_wall_model);
    const bool actual_reference = (m_reference_distance == "actual");
    const bool center_weight = (m_drag_weight == "center");

    WallParams wp{};
    wp.kappa = m_kappa;
    wp.beta_m = m_beta_m;
    wp.gamma_m = m_gamma_m;
    wp.inv_mo_length =
        (m_wall_het_model == "mol") ? 1.0_rt / m_monin_obukhov_length : 0.0_rt;
    wp.dt = dt;
    wp.tau_floor = m_forcing_time_factor * dt;
    wp.use_wall_time_scale = (m_bc_forcing_time_scale == "wall");

    amrex::ParallelFor(
        src_term, [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
            const auto& vel = vel_arrs[nbx];
            const auto& frac = frac_arrs[nbx];
            const auto& src = src_arrs[nbx];

            const amrex::Real beta = frac(i, j, k, 0);
            const amrex::Real w_solid = kynema_sgf::immersed_wall::solid_weight(
                beta, center_weight, solid_threshold);
            const int cell_mask = mask_arrs[nbx](i, j, k, 0);
            const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> u{
                vel(i, j, k, 0), vel(i, j, k, 1), vel(i, j, k, 2)};

            // 1. Immersed drag toward zero velocity, exact in time
            if (w_solid > 0.0_rt && drag_rate > 0.0_rt) {
                const amrex::Real C_eff =
                    kynema_sgf::immersed_wall::exact_relaxation_rate(
                        w_solid * drag_rate, dt);
                for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                    src(i, j, k, n) -= C_eff * u[n];
                }
            }

            // 2. Wall model in surface cells
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
                solid_threshold, actual_reference, patches.data());

            amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> force{
                0.0_rt, 0.0_rt, 0.0_rt};
            amrex::Real weight_sum = 0.0_rt;
            for (int ip = 0; ip < np; ++ip) {
                const WallPatch& p = patches[ip];
                const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> u_ref{
                    vel(p.ir, p.jr, p.kr, 0), vel(p.ir, p.jr, p.kr, 1),
                    vel(p.ir, p.jr, p.kr, 2)};
                const auto ut_ref =
                    kynema_sgf::immersed_wall::tangential(u_ref, p.nrm);
                const auto ut = kynema_sgf::immersed_wall::tangential(u, p.nrm);
                amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> f_patch{};
                kynema_sgf::immersed_wall::wall_momentum_forcing(
                    ut_ref, ut, p, z0, wp, f_patch);
                for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                    force[n] += p.weight * f_patch[n];
                }
                weight_sum += p.weight;
            }
            if (weight_sum > 0.0_rt) {
                for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                    src(i, j, k, n) +=
                        (1.0_rt - w_solid) * force[n] / weight_sum;
                }
            }
        });
    amrex::Gpu::streamSynchronize();
}

} // namespace kynema_sgf::pde::icns
