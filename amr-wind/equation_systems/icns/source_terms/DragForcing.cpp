#include "amr-wind/equation_systems/icns/source_terms/DragForcing.H"
#include "amr-wind/utilities/IOManager.H"

#include "AMReX_Gpu.H"
#include "AMReX_Random.H"
#include "amr-wind/wind_energy/ABL.H"

namespace amr_wind::pde::icns {

DragForcing::DragForcing(const CFDSim& sim)
    : m_sim(sim)
    , m_mesh(sim.mesh())
    , m_velocity(sim.repo().get_field("velocity"))
{
    amrex::ParmParse pp("DragForcing");
    pp.query("drag_coefficient", m_drag);
    pp.query("sponge_strength", m_sponge_strength);
    pp.query("sponge_density", m_sponge_density);
    pp.query("sponge_distanceX", m_sponge_distanceX);
    pp.query("sponge_distanceY", m_sponge_distanceY);
    const auto& phy_mgr = m_sim.physics_manager();
    if (phy_mgr.contains("ABL")) {
        const auto& abl = m_sim.physics_manager().get<amr_wind::ABL>();
        const VelPlaneAveraging& fa_velocity =
            abl.abl_statistics().vel_profile_coarse();
        device_vel_ht.resize(fa_velocity.line_centroids().size());
        device_vel_vals.resize(fa_velocity.line_average().size());
        amrex::Gpu::copy(
            amrex::Gpu::hostToDevice, fa_velocity.line_centroids().begin(),
            fa_velocity.line_centroids().end(), device_vel_ht.begin());
        amrex::Gpu::copy(
            amrex::Gpu::hostToDevice, fa_velocity.line_average().begin(),
            fa_velocity.line_average().end(), device_vel_vals.begin());
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
    const bool is_terrain = this->m_sim.repo().field_exists("terrainBlank");
    if (!is_terrain) {
        amrex::Abort("Need terrain blanking variable to use this source term");
    }
    auto* const m_terrain_blank = &this->m_sim.repo().get_field("terrainBlank");
    const auto& blank = (*m_terrain_blank)(lev).const_array(mfi);
    auto* const m_terrain_drag = &this->m_sim.repo().get_field("terrainDrag");
    const auto& drag = (*m_terrain_drag)(lev).const_array(mfi);
    auto* const m_terrainz0 = &this->m_sim.repo().get_field("terrainz0");
    const auto& terrainz0 = (*m_terrainz0)(lev).const_array(mfi);
    const auto& geom_vec = m_mesh.Geom();
    const auto& geom = geom_vec[lev];
    const auto& dx = geom.CellSize();
    const auto& prob_lo = geom.ProbLoArray();
    const auto& prob_hi = geom.ProbHiArray();
    const amrex::Real device_drag = m_drag;
    const amrex::Real device_sponge_strength = m_sponge_strength;
    const amrex::Real device_sponge_density = m_sponge_density;
    const amrex::Real device_startX = (m_sponge_distanceX > 0)
                                          ? prob_hi[0] - m_sponge_distanceX
                                          : prob_lo[0] - m_sponge_distanceX;
    const amrex::Real device_startY = (m_sponge_distanceY > 0)
                                          ? prob_hi[1] - m_sponge_distanceY
                                          : prob_lo[1] - m_sponge_distanceY;
    const unsigned device_spongeX = (m_sponge_distanceX > 0) ? 1 : 0;
    const unsigned device_spongeY = (m_sponge_distanceY > 0) ? 1 : 0;
    // Copy Data
    const auto* local_device_vel_ht = device_vel_ht.data();
    const auto* local_device_vel_vals = device_vel_vals.data();
    const unsigned vsize = device_vel_ht.size();
    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        const amrex::Real x1 = prob_lo[0] + (i + 0.5) * dx[0];
        const amrex::Real x2 = prob_lo[1] + (j + 0.5) * dx[1];
        const amrex::Real x3 = prob_lo[2] + (k + 0.5) * dx[2];
        amrex::Real xdamping = 0;
        amrex::Real ydamping = 0;
        amrex::Real xi =
            (device_spongeX == 1)
                ? (x1 - device_startX) / (prob_hi[0] - device_startX)
                : (device_startX - x1) / (device_startX - prob_lo[0]);
        xi = std::max(xi, 0.0);
        xdamping = device_sponge_strength * xi * xi;
        amrex::Real yi =
            (device_spongeY == 1)
                ? (x2 - device_startY) / (prob_hi[1] - device_startY)
                : (device_startY - x2) / (device_startY - prob_lo[1]);
        yi = std::max(yi, 0.0);
        ydamping = device_sponge_strength * yi * yi;
        const amrex::Real Cd = device_drag / dx[0];
        amrex::Real device_spongeVelX = 0.0;
        amrex::Real device_spongeVelY = 0.0;
        amrex::Real device_spongeVelZ = 0.0;
        amrex::Real residual = 1000;
        amrex::Real height_error = 0.0;
        for (unsigned ii = 0; ii < vsize; ++ii) {
            height_error = std::abs(x3 - local_device_vel_ht[ii]);
            if (height_error < residual) {
                residual = height_error;
                const unsigned ix = 3 * ii;
                const unsigned iy = 3 * ii + 1;
                const unsigned iz = 3 * ii + 2;
                device_spongeVelX = local_device_vel_vals[ix];
                device_spongeVelY = local_device_vel_vals[iy];
                device_spongeVelZ = local_device_vel_vals[iz];
            }
        }
        // Terrain Drag
        amrex::Real Dxz = 0.0;
        amrex::Real Dyz = 0.0;
        const amrex::Real ux1 = vel(i, j, k, 0);
        const amrex::Real uy1 = vel(i, j, k, 1);
        const amrex::Real uz1 = vel(i, j, k, 2);
        const amrex::Real m = std::sqrt(ux1 * ux1 + uy1 * uy1 + uz1 * uz1);
        if (drag(i, j, k) == 1) {
            const amrex::Real m1 = std::sqrt(ux1 * ux1 + uy1 * uy1);
            const amrex::Real ux2 = vel(i, j, k - 1, 0);
            const amrex::Real uy2 = vel(i, j, k - 1, 1);
            const amrex::Real m2 = std::sqrt(ux2 * ux2 + uy2 * uy2);
            const amrex::Real kappa = 0.41;
            const amrex::Real z0 = std::max(terrainz0(i, j, k), 1e-4);
            const amrex::Real ustar = std::min(device_sponge_strength, 1.0) *
                                      std::abs(m2 - m1) * kappa /
                                      std::log((x3 + z0) / z0);
            Dxz = -ustar * ustar * ux1 /
                  (1e-5 + std::sqrt(ux1 * ux1 + uy1 * uy1)) / dx[2];
            Dyz = -ustar * ustar * uy1 /
                  (1e-5 + std::sqrt(ux1 * ux1 + uy1 * uy1)) / dx[2];
        }
        // Adjusting Cd for momentum
        const amrex::Real CdM = std::min(Cd / (m + 1e-5), 1000.0);
        src_term(i, j, k, 0) -=
            (CdM * m * ux1 * blank(i, j, k) + Dxz * drag(i, j, k) +
             (xdamping + ydamping) *
                 (ux1 - device_sponge_density * device_spongeVelX));
        src_term(i, j, k, 1) -=
            (CdM * m * uy1 * blank(i, j, k) + Dyz * drag(i, j, k) +
             (xdamping + ydamping) *
                 (uy1 - device_sponge_density * device_spongeVelY));
        src_term(i, j, k, 2) -=
            (CdM * m * uz1 * blank(i, j, k) +
             (xdamping + ydamping) *
                 (uz1 - device_sponge_density * device_spongeVelZ));
    });
}

} // namespace amr_wind::pde::icns
