#include "src/equation_systems/icns/source_terms/ImmersedDragForcing.H"
#include "src/physics/ImmersedTerrain.H"
#include "src/utilities/constants.H"
#include "src/utilities/trig_ops.H"
#include "AMReX_Gpu.H"
#include "AMReX_ParmParse.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace {

using kynema_sgf::immersedterrain::ImmersedTerrain;

enum class WallModel : int { cell_offset = 0, terrain_height, surface_normal };

//! Monin-Obukhov stability function for momentum (device copy of
//! MOData::calc_psi_m, which is host only)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real psi_m(
    const amrex::Real zeta, const amrex::Real beta_m, const amrex::Real gamma_m)
{
    if (zeta > 0.0_rt) {
        return -gamma_m * zeta;
    }
    const amrex::Real x = std::sqrt(std::sqrt(1.0_rt - (beta_m * zeta)));
    return (2.0_rt * std::log(0.5_rt * (1.0_rt + x))) +
           std::log(0.5_rt * (1.0_rt + (x * x))) - (2.0_rt * std::atan(x)) +
           kynema_sgf::utils::half_pi();
}

//! Exact integration of du/dt = -C (u - u_t) over one step: the source
//! coefficient that reproduces u_t + (u - u_t) exp(-C dt)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
exact_relaxation_rate(const amrex::Real C, const amrex::Real dt)
{
    return (1.0_rt - std::exp(-C * dt)) / dt;
}

//! Constants shared by all wall-model variants
struct WallParams
{
    amrex::Real kappa;
    amrex::Real beta_m;
    amrex::Real gamma_m;
    amrex::Real inv_mo_length; // 1/L, zero for neutral
    amrex::Real dt;
    amrex::Real tau_floor; // tau_f * dt
    bool use_wall_time_scale;
};

/** Log-law forcing for one wall orientation.
 *
 *  \param ut_ref  tangential reference velocity at distance d2 from the wall
 *  \param ut      tangential velocity of the forced cell at distance d1
 *  \param d1, d2  wall distances of the cell and the reference point
 *  \param dxn     cell width along the wall normal
 *  \param z0      roughness
 *  \param[out] force  forcing on the tangential velocity (same frame as ut)
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void wall_forcing(
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& ut_ref,
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& ut,
    const amrex::Real d1,
    const amrex::Real d2,
    const amrex::Real dxn,
    const amrex::Real z0,
    const WallParams& wp,
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& force)
{
    const amrex::Real m_ref = std::sqrt(
        (ut_ref[0] * ut_ref[0]) + (ut_ref[1] * ut_ref[1]) +
        (ut_ref[2] * ut_ref[2]));
    const amrex::Real dd1 = amrex::max<amrex::Real>(d1, z0);
    const amrex::Real dd2 = amrex::max<amrex::Real>(d2, dd1);
    const amrex::Real phi2 =
        std::log(dd2 / z0) -
        psi_m(dd2 * wp.inv_mo_length, wp.beta_m, wp.gamma_m);
    const amrex::Real phi1 =
        std::log(dd1 / z0) -
        psi_m(dd1 * wp.inv_mo_length, wp.beta_m, wp.gamma_m);
    const amrex::Real ustar = wp.kappa * m_ref / phi2;
    const amrex::Real ut1_mag = ustar / wp.kappa * phi1;

    const amrex::Real tau_wall =
        wp.use_wall_time_scale
            ? amrex::max<amrex::Real>(
                  wp.tau_floor, dd1 / (ustar + kynema_sgf::constants::EPS))
            : wp.tau_floor;
    const amrex::Real C_bc = exact_relaxation_rate(1.0_rt / tau_wall, wp.dt);

    for (int n = 0; n < AMREX_SPACEDIM; ++n) {
        const amrex::Real e_t =
            ut_ref[n] / (m_ref + kynema_sgf::constants::EPS);
        const amrex::Real stress_div = -ustar * ustar * e_t / dxn;
        const amrex::Real relax = -C_bc * (ut[n] - (ut1_mag * e_t));
        force[n] = stress_div + relax;
    }
}

} // namespace

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
    pp.query("minimum_z0", m_min_z0);
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

    const amrex::Real drag_rate = m_drag_coefficient / dx[2];
    const amrex::Real min_z0 = m_min_z0;
    const amrex::Real solid_threshold = m_solid_threshold;
    const bool apply_wall_model = !m_is_laminar;
    const WallModel wall_model =
        (m_wall_model == "surface_normal")   ? WallModel::surface_normal
        : (m_wall_model == "terrain_height") ? WallModel::terrain_height
                                             : WallModel::cell_offset;

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
            const auto& surf = surf_arrs[nbx];
            auto& src = src_arrs[nbx];

            const amrex::Real beta = frac(i, j, k, 0);
            const int cell_mask = mask_arrs[nbx](i, j, k, 0);
            amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> u{
                vel(i, j, k, 0), vel(i, j, k, 1), vel(i, j, k, 2)};

            // 1. Immersed drag toward zero velocity, exact in time
            if (beta > 0.0_rt) {
                const amrex::Real C_eff =
                    exact_relaxation_rate(beta * drag_rate, dt);
                for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                    src(i, j, k, n) -= C_eff * u[n];
                }
            }

            // 2. Wall model in surface cells
            if (!apply_wall_model ||
                cell_mask != ImmersedTerrain::mask_surface || beta >= 1.0_rt) {
                return;
            }
            const amrex::Real z0 =
                amrex::max<amrex::Real>(z0_arrs[nbx](i, j, k, 0), min_z0);
            const amrex::Real z_c = prob_lo[2] + ((k + 0.5_rt) * dx[2]);
            const amrex::Real h = surf(i, j, k, ImmersedTerrain::surf_height);

            amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> force{
                0.0_rt, 0.0_rt, 0.0_rt};

            // Terrain below the cell: partial cell or solid bottom neighbour
            const bool terrain_below =
                (beta > 0.0_rt) || (frac(i, j, k - 1, 0) >= solid_threshold);

            if (wall_model == WallModel::surface_normal && terrain_below) {
                // Unit normal from the surface slopes, pointing into the air
                const amrex::Real hx =
                    surf(i, j, k, ImmersedTerrain::surf_slope_x);
                const amrex::Real hy =
                    surf(i, j, k, ImmersedTerrain::surf_slope_y);
                const amrex::Real nmag =
                    std::sqrt(1.0_rt + (hx * hx) + (hy * hy));
                const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> nrm{
                    -hx / nmag, -hy / nmag, 1.0_rt / nmag};
                // Cell width along the normal and wall distances
                const amrex::Real dxn = 1.0_rt / ((std::abs(nrm[0]) / dx[0]) +
                                                  (std::abs(nrm[1]) / dx[1]) +
                                                  (std::abs(nrm[2]) / dx[2]));
                const amrex::Real d1 = (z_c - h) * nrm[2];
                const amrex::Real d2 = amrex::max<amrex::Real>(d1, z0) + dxn;
                // Nearest cell one normal step away, within the ghost layer
                const auto offset = [=](const int dir) {
                    const amrex::Real s = dxn * nrm[dir] / dx[dir];
                    const int o = static_cast<int>(std::floor(s + 0.5_rt));
                    return amrex::min(1, amrex::max(-1, o));
                };
                const int ir = i + offset(0);
                const int jr = j + offset(1);
                const int kr = k + offset(2);
                amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> u_ref{
                    vel(ir, jr, kr, 0), vel(ir, jr, kr, 1), vel(ir, jr, kr, 2)};
                // Tangential projections
                const amrex::Real un_ref = (u_ref[0] * nrm[0]) +
                                           (u_ref[1] * nrm[1]) +
                                           (u_ref[2] * nrm[2]);
                const amrex::Real un =
                    (u[0] * nrm[0]) + (u[1] * nrm[1]) + (u[2] * nrm[2]);
                amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> ut_ref{};
                amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> ut{};
                for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                    ut_ref[n] = u_ref[n] - (un_ref * nrm[n]);
                    ut[n] = u[n] - (un * nrm[n]);
                }
                wall_forcing(ut_ref, ut, d1, d2, dxn, z0, wp, force);
            } else {
                // Face based: loop over the six faces, accumulate the forcing
                // from every face that touches a solid neighbour, weighted by
                // that neighbour's terrain fraction
                amrex::Real weight_sum = 0.0_rt;
                for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
                    for (int s = -1; s <= 1; s += 2) {
                        const int in = i + ((dir == 0) ? s : 0);
                        const int jn = j + ((dir == 1) ? s : 0);
                        const int kn = k + ((dir == 2) ? s : 0);
                        const amrex::Real beta_nb = frac(in, jn, kn, 0);
                        if (beta_nb < solid_threshold) {
                            continue;
                        }
                        // Reference cell is the neighbour opposite the wall
                        const int ir = i - ((dir == 0) ? s : 0);
                        const int jr = j - ((dir == 1) ? s : 0);
                        const int kr = k - ((dir == 2) ? s : 0);

                        amrex::Real d1 = 0.5_rt * dx[dir];
                        amrex::Real d2 = 1.5_rt * dx[dir];
                        if (wall_model == WallModel::terrain_height &&
                            dir == 2 && s == -1) {
                            d1 = z_c - h;
                            d2 = z_c + dx[2] - h;
                        }
                        // Tangential components are the two not along dir
                        amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> ut_ref{};
                        amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> ut{};
                        for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                            ut_ref[n] =
                                (n == dir) ? 0.0_rt : vel(ir, jr, kr, n);
                            ut[n] = (n == dir) ? 0.0_rt : u[n];
                        }
                        amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> f_face{};
                        wall_forcing(
                            ut_ref, ut, d1, d2, dx[dir], z0, wp, f_face);
                        for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                            force[n] += beta_nb * f_face[n];
                        }
                        weight_sum += beta_nb;
                    }
                }
                if (weight_sum > 0.0_rt) {
                    for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                        force[n] /= weight_sum;
                    }
                }
            }

            for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                src(i, j, k, n) += (1.0_rt - beta) * force[n];
            }
        });
    amrex::Gpu::streamSynchronize();
}

} // namespace kynema_sgf::pde::icns
