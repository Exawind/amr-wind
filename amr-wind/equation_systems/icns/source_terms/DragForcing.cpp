#include "amr-wind/equation_systems/icns/source_terms/DragForcing.H"
#include "amr-wind/utilities/IOManager.H"

#include "AMReX_Gpu.H"
#include "AMReX_Random.H"
#include "amr-wind/wind_energy/ABL.H"
#include "amr-wind/utilities/linear_interpolation.H"

namespace amr_wind::pde::icns {

DragForcing::DragForcing(const CFDSim& sim)
    : m_time(sim.time())
    , m_sim(sim)
    , m_mesh(sim.mesh())
    , m_velocity(sim.repo().get_field("velocity"))
{
    amrex::ParmParse pp("DragForcing");
    pp.query("drag_coefficient", m_drag_coefficient);
    pp.query("sponge_strength", m_sponge_strength);
    pp.query("sponge_density", m_sponge_density);
    pp.query("sponge_distance_west", m_sponge_distance_west);
    pp.query("sponge_distance_east", m_sponge_distance_east);
    pp.query("sponge_distance_south", m_sponge_distance_south);
    pp.query("sponge_distance_north", m_sponge_distance_north);
    pp.query("sponge_west", m_sponge_west);
    pp.query("sponge_east", m_sponge_east);
    pp.query("sponge_south", m_sponge_south);
    pp.query("sponge_north", m_sponge_north);
    pp.query("is_laminar", m_is_laminar);
    const auto& phy_mgr = m_sim.physics_manager();
    if (phy_mgr.contains("ABL")) {
        const auto& abl = m_sim.physics_manager().get<amr_wind::ABL>();
        const auto& fa_velocity = abl.abl_statistics().vel_profile_coarse();
        m_device_vel_ht.resize(fa_velocity.line_centroids().size());
        m_device_vel_vals.resize(fa_velocity.line_average().size());
        amrex::Gpu::copy(
            amrex::Gpu::hostToDevice, fa_velocity.line_centroids().begin(),
            fa_velocity.line_centroids().end(), m_device_vel_ht.begin());
        amrex::Gpu::copy(
            amrex::Gpu::hostToDevice, fa_velocity.line_average().begin(),
            fa_velocity.line_average().end(), m_device_vel_vals.begin());
    } else {
        m_sponge_strength = 0.0;
    }
}

DragForcing::~DragForcing() = default;

void DragForcing::operator()(
    const int lev,
    const amrex::MFIter& mfi,
    const amrex::Box& bx,
    const FieldState fstate,
    const amrex::Array4<amrex::Real>& src_term) const
{
    const auto& vel =
        m_velocity.state(field_impl::dof_state(fstate))(lev).const_array(mfi);
    const bool is_terrain =
        this->m_sim.repo().int_field_exists("terrain_blank");
    if (!is_terrain) {
        amrex::Abort("Need terrain blanking variable to use this source term");
    }
    auto* const m_terrain_drag =
        &this->m_sim.repo().get_int_field("terrain_drag");
    const auto& drag = (*m_terrain_drag)(lev).const_array(mfi);
    auto* const m_terrainz0 = &this->m_sim.repo().get_field("terrainz0");
    const auto& terrainz0 = (*m_terrainz0)(lev).const_array(mfi);
    const auto* m_terrain_vf = &this->m_sim.repo().get_field("terrain_vf");
    const auto& vf_arr = (*m_terrain_vf)(lev).const_array(mfi);
    const auto& geom = m_mesh.Geom(lev);
    const auto& dx = geom.CellSizeArray();
    const auto& prob_lo = geom.ProbLoArray();
    const auto& prob_hi = geom.ProbHiArray();
    const amrex::Real drag_coefficient = m_drag_coefficient;
    const amrex::Real sponge_strength = m_sponge_strength;
    const amrex::Real sponge_density = m_sponge_density;
    const amrex::Real start_east = prob_hi[0] - m_sponge_distance_east;
    const amrex::Real start_west = prob_lo[0] - m_sponge_distance_west;
    const amrex::Real start_north = prob_hi[1] - m_sponge_distance_north;
    const amrex::Real start_south = prob_lo[1] - m_sponge_distance_south;
    const int sponge_east = m_sponge_east;
    const int sponge_west = m_sponge_west;
    const int sponge_south = m_sponge_south;
    const int sponge_north = m_sponge_north;
    // Copy Data
    const auto* device_vel_ht = m_device_vel_ht.data();
    const auto* device_vel_vals = m_device_vel_vals.data();
    const unsigned vsize = m_device_vel_ht.size();
    const auto& dt = m_time.delta_t();
    const bool is_laminar = m_is_laminar;
    const amrex::Real scale_factor = (dx[2] < 1.0) ? 1.0 : 1.0 / dx[2];
    const amrex::Real Cd =
        (is_laminar && dx[2] < 1) ? drag_coefficient : drag_coefficient / dx[2];
    const amrex::Real z0_min = 1e-4;
    const auto tiny = std::numeric_limits<amrex::Real>::epsilon();
    const amrex::Real kappa = 0.41;
    const amrex::Real cd_max = 1000.0;
    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        const amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];
        const amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1];
        const amrex::Real z = prob_lo[2] + (k + 0.5) * dx[2];
        amrex::Real xstart_damping = 0;
        amrex::Real ystart_damping = 0;
        amrex::Real xend_damping = 0;
        amrex::Real yend_damping = 0;
        amrex::Real xi_end = (x - start_east) / (prob_hi[0] - start_east);
        amrex::Real xi_start = (start_west - x) / (start_west - prob_lo[0]);
        xi_start = sponge_west * std::max(xi_start, 0.0);
        xi_end = sponge_east * std::max(xi_end, 0.0);
        xstart_damping = sponge_west * sponge_strength * xi_start * xi_start;
        xend_damping = sponge_east * sponge_strength * xi_end * xi_end;
        amrex::Real yi_end = (y - start_north) / (prob_hi[1] - start_north);
        amrex::Real yi_start = (start_south - y) / (start_south - prob_lo[1]);
        yi_start = sponge_south * std::max(yi_start, 0.0);
        yi_end = sponge_north * std::max(yi_end, 0.0);
        ystart_damping = sponge_strength * yi_start * yi_start;
        yend_damping = sponge_strength * yi_end * yi_end;
        const amrex::Real ux1 = vel(i, j, k, 0);
        const amrex::Real uy1 = vel(i, j, k, 1);
        const amrex::Real uz1 = vel(i, j, k, 2);
        const auto idx =
            interp::bisection_search(device_vel_ht, device_vel_ht + vsize, z);
        const amrex::Real spongeVelX =
            (vsize > 0) ? interp::linear_impl(
                              device_vel_ht, device_vel_vals, z, idx, 3, 0)
                        : ux1;
        const amrex::Real spongeVelY =
            (vsize > 0) ? interp::linear_impl(
                              device_vel_ht, device_vel_vals, z, idx, 3, 1)
                        : uy1;
        const amrex::Real spongeVelZ =
            (vsize > 0) ? interp::linear_impl(
                              device_vel_ht, device_vel_vals, z, idx, 3, 2)
                        : uz1;
        amrex::Real Dxz = 0.0;
        amrex::Real Dyz = 0.0;
        amrex::Real bc_forcing_x = 0;
        amrex::Real bc_forcing_y = 0;

        const amrex::Real m = std::sqrt(ux1 * ux1 + uy1 * uy1 + uz1 * uz1);
        if (drag(i, j, k) == 1 && (!is_laminar)) {
            const amrex::Real ux2 = vel(i, j, k + 1, 0);
            const amrex::Real uy2 = vel(i, j, k + 1, 1);
            const amrex::Real m2 = std::sqrt(ux2 * ux2 + uy2 * uy2);
            const amrex::Real z0 = std::max(terrainz0(i, j, k), z0_min);
            const amrex::Real ustar = m2 * kappa / std::log(1.5 * dx[2] / z0);
            const amrex::Real uTarget =
                ustar / kappa * std::log(0.5 * dx[2] / z0);
            const amrex::Real uxTarget =
                uTarget * ux2 / (tiny + std::sqrt(ux2 * ux2 + uy2 * uy2));
            const amrex::Real uyTarget =
                uTarget * uy2 / (tiny + std::sqrt(ux2 * ux2 + uy2 * uy2));
            bc_forcing_x = -(uxTarget - ux1) / dt;
            bc_forcing_y = -(uyTarget - uy1) / dt;
            Dxz = -ustar * ustar * ux1 /
                  (tiny + std::sqrt(ux1 * ux1 + uy1 * uy1)) / dx[2];
            Dyz = -ustar * ustar * uy1 /
                  (tiny + std::sqrt(ux1 * ux1 + uy1 * uy1)) / dx[2];
        }
        const amrex::Real CdM =
            std::min(Cd / (m + tiny), cd_max / scale_factor);
        src_term(i, j, k, 0) -=
            (CdM * m * ux1 * vf_arr(i, j, k) + Dxz * drag(i, j, k) +
             bc_forcing_x * drag(i, j, k) +
             (xstart_damping + xend_damping + ystart_damping + yend_damping) *
                 (ux1 - sponge_density * spongeVelX));
        src_term(i, j, k, 1) -=
            (CdM * m * uy1 * vf_arr(i, j, k) + Dyz * drag(i, j, k) +
             bc_forcing_y * drag(i, j, k) +
             (xstart_damping + xend_damping + ystart_damping + yend_damping) *
                 (uy1 - sponge_density * spongeVelY));
        src_term(i, j, k, 2) -=
            (CdM * m * uz1 * vf_arr(i, j, k) +
             (xstart_damping + xend_damping + ystart_damping + yend_damping) *
                 (uz1 - sponge_density * spongeVelZ));
    });
}

} // namespace amr_wind::pde::icns
