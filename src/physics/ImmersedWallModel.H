#ifndef IMMERSEDWALLMODEL_H
#define IMMERSEDWALLMODEL_H

#include "src/physics/ImmersedTerrain.H"
#include "src/utilities/constants.H"
#include "src/utilities/trig_ops.H"
#include "AMReX_Array4.H"
#include "AMReX_GpuQualifiers.H"
#include "AMReX_REAL.H"

#include <cstdint>
#include <string>

using namespace amrex::literals;

/** Wall-model geometry and Monin-Obukhov helpers shared by the immersed
 *  terrain source terms (ImmersedDragForcing, ImmersedDragTempForcing).
 *
 *  A surface cell (terrain_mask == 2) receives one or more "wall patches".
 *  Each patch carries the wall distance of the cell center d1, the distance
 *  d2 of a reference cell whose velocity and temperature are trusted, the cell
 *  width along the wall normal dxn, the unit normal pointing into the fluid,
 *  the reference cell indices and a weight. Both source terms build the same
 *  patches, so momentum and heat see the same friction velocity.
 */
namespace kynema_sgf::immersed_wall {

enum class WallModel : std::uint8_t {
    cell_offset = 0,
    terrain_height,
    surface_normal
};

inline WallModel parse_wall_model(const std::string& name)
{
    if (name == "surface_normal") {
        return WallModel::surface_normal;
    }
    if (name == "terrain_height") {
        return WallModel::terrain_height;
    }
    return WallModel::cell_offset;
}

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

//! Monin-Obukhov stability function for heat (device copy of
//! MOData::calc_psi_h)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real psi_h(
    const amrex::Real zeta, const amrex::Real beta_h, const amrex::Real gamma_h)
{
    if (zeta > 0.0_rt) {
        return -gamma_h * zeta;
    }
    const amrex::Real x = std::sqrt(1.0_rt - (beta_h * zeta));
    return 2.0_rt * std::log(0.5_rt * (1.0_rt + x));
}

//! Exact integration of du/dt = -C (u - u_t) over one step: the source
//! coefficient that reproduces u_t + (u - u_t) exp(-C dt)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
exact_relaxation_rate(const amrex::Real C, const amrex::Real dt)
{
    return (1.0_rt - std::exp(-C * dt)) / dt;
}

//! Constants shared by the wall-model variants
struct WallParams
{
    amrex::Real kappa{0.41_rt};
    amrex::Real beta_m{16.0_rt};
    amrex::Real gamma_m{5.0_rt};
    amrex::Real beta_h{16.0_rt};
    amrex::Real gamma_h{5.0_rt};
    amrex::Real inv_mo_length{0.0_rt}; //!< 1/L, zero for neutral
    amrex::Real dt{0.0_rt};
    amrex::Real tau_floor{0.0_rt}; //!< tau_f * dt
    bool use_wall_time_scale{true};

    //! ln(d/z0) - psi_m(d/L)
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
    phi_m(const amrex::Real d, const amrex::Real z0) const
    {
        return std::log(d / z0) - psi_m(d * inv_mo_length, beta_m, gamma_m);
    }
    //! ln(d/z0) - psi_h(d/L)
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
    phi_h(const amrex::Real d, const amrex::Real z0) const
    {
        return std::log(d / z0) - psi_h(d * inv_mo_length, beta_h, gamma_h);
    }
    //! Relaxation rate toward the wall-model target for a cell at distance d1
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
    relaxation_rate(const amrex::Real d1, const amrex::Real ustar) const
    {
        const amrex::Real tau =
            use_wall_time_scale
                ? amrex::max<amrex::Real>(
                      tau_floor, d1 / (ustar + kynema_sgf::constants::EPS))
                : tau_floor;
        return exact_relaxation_rate(1.0_rt / tau, dt);
    }
};

//! One wall contribution to a surface cell
struct WallPatch
{
    amrex::Real d1{0.0_rt};     //!< wall distance of the cell center
    amrex::Real d2{0.0_rt};     //!< wall distance of the reference cell
    amrex::Real dxn{0.0_rt};    //!< cell width along the wall normal
    amrex::Real weight{0.0_rt}; //!< weight when several patches are active
    int ir{0}, jr{0}, kr{0};    //!< reference cell
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> nrm{0.0_rt, 0.0_rt, 1.0_rt};
};

//! Component of v tangential to the patch normal
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> tangential(
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& v,
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& nrm)
{
    const amrex::Real vn = (v[0] * nrm[0]) + (v[1] * nrm[1]) + (v[2] * nrm[2]);
    return {v[0] - (vn * nrm[0]), v[1] - (vn * nrm[1]), v[2] - (vn * nrm[2])};
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
magnitude(const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& v)
{
    return std::sqrt((v[0] * v[0]) + (v[1] * v[1]) + (v[2] * v[2]));
}

/** Build the wall patches of a surface cell.
 *
 *  Face based (cell_offset, terrain_height): one patch per face whose
 *  neighbor has a terrain fraction at or above solid_threshold. The
 *  reference cell is the neighbor opposite the wall, d1 = 0.5 dx_f,
 *  d2 = 1.5 dx_f, except that terrain_height uses z - h on the bottom face.
 *  Patches are weighted by the neighbor fraction.
 *
 *  Surface normal: a single patch along the terrain normal built from the
 *  slopes, with the reference cell one normal cell-width away (nearest cell,
 *  offset clamped to the one available ghost layer). Cells with no terrain
 *  below (side walls) fall back to the face-based patches.
 *
 *  \return number of patches written (at most 6)
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE int wall_patches(
    const WallModel model,
    const int i,
    const int j,
    const int k,
    const amrex::Real beta,
    amrex::Array4<amrex::Real const> const& frac,
    amrex::Array4<amrex::Real const> const& surf,
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& dx,
    const amrex::Real z_c,
    const amrex::Real z0,
    const amrex::Real solid_threshold,
    WallPatch* patches)
{
    using kynema_sgf::immersedterrain::ImmersedTerrain;
    const amrex::Real h = surf(i, j, k, ImmersedTerrain::surf_height);
    const bool terrain_below =
        (beta > 0.0_rt) || (frac(i, j, k - 1, 0) >= solid_threshold);

    if (model == WallModel::surface_normal && terrain_below) {
        const amrex::Real hx = surf(i, j, k, ImmersedTerrain::surf_slope_x);
        const amrex::Real hy = surf(i, j, k, ImmersedTerrain::surf_slope_y);
        const amrex::Real nmag = std::sqrt(1.0_rt + (hx * hx) + (hy * hy));
        WallPatch& p = patches[0];
        p.nrm = {-hx / nmag, -hy / nmag, 1.0_rt / nmag};
        p.dxn = 1.0_rt /
                ((std::abs(p.nrm[0]) / dx[0]) + (std::abs(p.nrm[1]) / dx[1]) +
                 (std::abs(p.nrm[2]) / dx[2]));
        p.d1 = (z_c - h) * p.nrm[2];
        p.d2 = amrex::max<amrex::Real>(p.d1, z0) + p.dxn;
        const auto offset = [&](const int dir) {
            const amrex::Real s = p.dxn * p.nrm[dir] / dx[dir];
            const int o = static_cast<int>(std::floor(s + 0.5_rt));
            return amrex::min(1, amrex::max(-1, o));
        };
        p.ir = i + offset(0);
        p.jr = j + offset(1);
        p.kr = k + offset(2);
        p.weight = 1.0_rt;
        return 1;
    }

    int np = 0;
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        for (int s = -1; s <= 1; s += 2) {
            const int in = i + ((dir == 0) ? s : 0);
            const int jn = j + ((dir == 1) ? s : 0);
            const int kn = k + ((dir == 2) ? s : 0);
            const amrex::Real beta_nb = frac(in, jn, kn, 0);
            if (beta_nb < solid_threshold) {
                continue;
            }
            WallPatch& p = patches[np++];
            p.ir = i - ((dir == 0) ? s : 0);
            p.jr = j - ((dir == 1) ? s : 0);
            p.kr = k - ((dir == 2) ? s : 0);
            // Normal points from the wall into the cell
            p.nrm = {0.0_rt, 0.0_rt, 0.0_rt};
            p.nrm[dir] = -static_cast<amrex::Real>(s);
            p.dxn = dx[dir];
            p.d1 = 0.5_rt * dx[dir];
            p.d2 = 1.5_rt * dx[dir];
            if (model == WallModel::terrain_height && dir == 2 && s == -1) {
                p.d1 = z_c - h;
                p.d2 = z_c + dx[2] - h;
            }
            p.weight = beta_nb;
        }
    }
    return np;
}

/** Friction velocity of a patch from the tangential reference speed */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real friction_velocity(
    const amrex::Real ut_ref_mag,
    const WallPatch& p,
    const amrex::Real z0,
    const WallParams& wp)
{
    const amrex::Real dd1 = amrex::max<amrex::Real>(p.d1, z0);
    const amrex::Real dd2 = amrex::max<amrex::Real>(p.d2, dd1);
    return wp.kappa * ut_ref_mag / wp.phi_m(dd2, z0);
}

/** Log-law momentum forcing of one patch: stress divergence plus relaxation
 *  of the tangential velocity toward the log-law value at d1.
 *
 *  \param ut_ref tangential reference velocity at d2
 *  \param ut     tangential velocity of the forced cell at d1
 *  \param p      wall patch (distances, normal, cell width)
 *  \param z0     roughness length
 *  \param wp     wall-model constants
 *  \param[out] force forcing on the tangential velocity (Cartesian frame)
 *  \return friction velocity
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real wall_momentum_forcing(
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& ut_ref,
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& ut,
    const WallPatch& p,
    const amrex::Real z0,
    const WallParams& wp,
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& force)
{
    const amrex::Real m_ref = magnitude(ut_ref);
    const amrex::Real dd1 = amrex::max<amrex::Real>(p.d1, z0);
    const amrex::Real ustar = friction_velocity(m_ref, p, z0, wp);
    const amrex::Real ut1_mag = ustar / wp.kappa * wp.phi_m(dd1, z0);
    const amrex::Real C_bc = wp.relaxation_rate(dd1, ustar);
    for (int n = 0; n < AMREX_SPACEDIM; ++n) {
        const amrex::Real e_t =
            ut_ref[n] / (m_ref + kynema_sgf::constants::EPS);
        force[n] =
            (-ustar * ustar * e_t / p.dxn) - (C_bc * (ut[n] - (ut1_mag * e_t)));
    }
    return ustar;
}

/** Log-law heat forcing of one patch given the temperature scale theta*:
 *  surface flux divergence plus relaxation toward the log-law temperature.
 *
 *  \param theta      temperature of the forced cell at d1
 *  \param theta_s    surface temperature
 *  \param theta_star temperature scale (positive for a cooling surface)
 *  \param ustar      friction velocity of the patch
 *  \param p          wall patch (distances, normal, cell width)
 *  \param z0         roughness length
 *  \param wp         wall-model constants
 *  \return forcing on the temperature
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real wall_heat_forcing(
    const amrex::Real theta,
    const amrex::Real theta_s,
    const amrex::Real theta_star,
    const amrex::Real ustar,
    const WallPatch& p,
    const amrex::Real z0,
    const WallParams& wp)
{
    const amrex::Real dd1 = amrex::max<amrex::Real>(p.d1, z0);
    const amrex::Real theta_target =
        theta_s + (theta_star / wp.kappa * wp.phi_h(dd1, z0));
    const amrex::Real q_s = -ustar * theta_star; // kinematic heat flux
    const amrex::Real C_bc = wp.relaxation_rate(dd1, ustar);
    return (q_s / p.dxn) - (C_bc * (theta - theta_target));
}

} // namespace kynema_sgf::immersed_wall

#endif
