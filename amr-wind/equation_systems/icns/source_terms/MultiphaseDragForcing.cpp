#include "amr-wind/equation_systems/icns/source_terms/MultiphaseDragForcing.H"
#include "amr-wind/equation_systems/vof/volume_fractions.H"
#include "amr-wind/utilities/IOManager.H"

#include "AMReX_Gpu.H"
#include "AMReX_Random.H"
#include "amr-wind/wind_energy/ABL.H"
#include "amr-wind/physics/TerrainDrag.H"
#include "amr-wind/utilities/linear_interpolation.H"

namespace {
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real viscous_drag_calculations(
    amrex::Real& Dxz,
    amrex::Real& Dyz,
    const amrex::Real ux1r,
    const amrex::Real uy1r,
    const amrex::Real ux2r,
    const amrex::Real uy2r,
    const amrex::Real z0,
    const amrex::Real dz,
    const amrex::Real kappa,
    const amrex::Real tiny)
{
    const amrex::Real m2 = std::sqrt(ux2r * ux2r + uy2r * uy2r);
    const amrex::Real ustar = m2 * kappa / std::log(1.5 * dz / z0);
    Dxz += -ustar * ustar * ux1r /
           (tiny + std::sqrt(ux1r * ux1r + uy1r * uy1r)) / dz;
    Dyz += -ustar * ustar * uy1r /
           (tiny + std::sqrt(ux1r * ux1r + uy1r * uy1r)) / dz;
    return ustar;
}

} // namespace

namespace amr_wind::pde::icns {

MultiphaseDragForcing::MultiphaseDragForcing(const CFDSim& sim)
    : m_time(sim.time())
    , m_sim(sim)
    , m_mesh(sim.mesh())
    , m_velocity(sim.repo().get_field("velocity"))
{
    const auto terrain_phys =
        m_sim.physics_manager().get<amr_wind::terraindrag::TerrainDrag>();
    const auto target_vel_name = terrain_phys.wave_velocity_field_name();
    m_target_vel = &sim.repo().get_field(target_vel_name);
    const auto target_levelset_name =
        terrain_phys.wave_negative_elevation_name();
    m_target_levelset = &sim.repo().get_field(target_levelset_name);
    const bool has_vof = sim.repo().field_exists("vof");
    if (!has_vof) {
        amrex::Abort("Works only when vof field is available");
    }
}

MultiphaseDragForcing::~MultiphaseDragForcing() = default;

void MultiphaseDragForcing::operator()(
    const int lev,
    const amrex::MFIter& mfi,
    const amrex::Box& bx,
    const FieldState fstate,
    const amrex::Array4<amrex::Real>& src_term) const
{
    const auto& vel =
        m_velocity.state(field_impl::dof_state(fstate))(lev).const_array(mfi);
    auto* const m_vof = &this->m_sim.repo().get_field("vof");
    const auto& vof = (*m_vof)(lev).const_array(mfi);
    auto* const m_terrainz0 = &this->m_sim.repo().get_field("terrainz0");
    const auto& terrainz0 = (*m_terrainz0)(lev).const_array(mfi);
    const auto& target_vel_arr = (*m_target_vel)(lev).const_array(mfi);
    const auto& target_lvs_arr = (*m_target_levelset)(lev).const_array(mfi);
    const auto& geom = m_mesh.Geom(lev);
    const auto& dx = geom.CellSizeArray();
    const auto& dt = m_time.delta_t();
    const amrex::Real kappa = 0.41;
    const amrex::Real z0_min = 1e-4;
    const auto tiny = std::numeric_limits<amrex::Real>::epsilon();
    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        const amrex::Real ux1 = vel(i, j, k, 0);
        const amrex::Real uy1 = vel(i, j, k, 1);
        amrex::Real Dxz = 0.0;
        amrex::Real Dyz = 0.0;
        amrex::Real bc_forcing_x = 0;
        amrex::Real bc_forcing_y = 0;
        int k_off = -1;
        if (vof(i, j, k) == 0 && vof(i, j, k - 1) > 0) {
            const amrex::Real cell_length_2D =
                std::sqrt(dx[0] * dx[0] + dx[2] * dx[2]);
            if (target_lvs_arr(i, j, k) + cell_length_2D >= 0) {
                // Current cell will be used for wave velocity
                k_off = 0;
            }
            const amrex::Real wall_u = target_vel_arr(i, j, k + k_off, 0);
            const amrex::Real wall_v = target_vel_arr(i, j, k + k_off, 1);
            // Relative velocities for calculating shear
            const amrex::Real ux1r = ux1 - wall_u;
            const amrex::Real uy1r = uy1 - wall_v;
            const amrex::Real ux2r = vel(i, j, k + 1, 0) - wall_u;
            const amrex::Real uy2r = vel(i, j, k + 1, 1) - wall_v;
            const amrex::Real z0 = std::max(terrainz0(i, j, k), z0_min);
            const amrex::Real ustar = viscous_drag_calculations(
                Dxz, Dyz, ux1r, uy1r, ux2r, uy2r, z0, dx[2], kappa, tiny);
            const amrex::Real uTarget =
                ustar / kappa * std::log(0.5 * dx[2] / z0);
            const amrex::Real uxTarget =
                uTarget * ux2r / (tiny + std::sqrt(ux2r * ux2r + uy2r * uy2r));
            const amrex::Real uyTarget =
                uTarget * uy2r / (tiny + std::sqrt(ux2r * ux2r + uy2r * uy2r));
            // BC forcing pushes nonrelative velocity toward target velocity
            bc_forcing_x = -(uxTarget - ux1) / dt;
            bc_forcing_y = -(uyTarget - uy1) / dt;
        }
        src_term(i, j, k, 0) -= Dxz + bc_forcing_x;
        src_term(i, j, k, 1) -= Dyz + bc_forcing_y;
    });
}

} // namespace amr_wind::pde::icns
