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
    , m_terrainBlank(sim.repo().get_field("terrainBlank"))
    , m_terrainDrag(sim.repo().get_field("terrainDrag"))

{
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
    const auto blank = m_terrainBlank(lev).const_array(mfi);
    const auto drag = m_terrainDrag(lev).const_array(mfi);
    const auto& geom_vec = m_mesh.Geom();
    const auto& geom = geom_vec[lev];
    const auto& dx = geom.CellSizeArray();
    const auto& prob_lo = geom.ProbLoArray();
    const auto& prob_hi = geom.ProbHiArray();
    const amrex::Real dragCoeff = m_dragCoeff;
    const amrex::Real spongeStrength = m_spongeStrength;
    const amrex::Real spongeDensity = m_spongeDensity;
    const amrex::Real startX = (1 - m_spongePercentX / 100.0) * prob_hi[0];
    const amrex::Real startY = (1 - m_spongePercentY / 100.0) * prob_hi[1];
    unsigned long verticalSize = device_vel_ht.size();
    amrex::Gpu::DeviceVector<amrex::Real> vel_ht;
    amrex::Gpu::DeviceVector<amrex::Real> vel_vals;
    vel_ht.resize(verticalSize);
    vel_vals.resize(verticalSize);
    amrex::Gpu::copy(
        amrex::Gpu::deviceToDevice, device_vel_ht.begin(),
        device_vel_ht.end(), vel_ht.begin());
    amrex::Gpu::copy(
        amrex::Gpu::deviceToDevice, device_vel_vals.begin(),
        device_vel_vals.end(), vel_vals.begin());
    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        const amrex::Real ux = vel(i, j, k, 0);
        const amrex::Real uy = vel(i, j, k, 1);
        const amrex::Real uz = vel(i, j, k, 2);
        const amrex::Real x1 = prob_lo[0] + (i + 0.5) * dx[0];
        const amrex::Real x2 = prob_lo[1] + (j + 0.5) * dx[1];
        const amrex::Real x3 = prob_lo[2] + (k + 0.5) * dx[2];
        amrex::Real xdamping = 0;
        amrex::Real ydamping = 0;
        if (x1 > startX) {
            amrex::Real xi = (x1 - startX) / (prob_hi[0] - startX);
            xdamping = spongeStrength * xi * xi;
        }
        if (x2 > startY) {
            amrex::Real yi = (x2 - startY) / (prob_hi[1] - startY);
            ydamping = spongeStrength * yi * yi;
        }
        const amrex::Real m = std::sqrt(ux * ux + uy * uy + uz * uz);
        amrex::Real Cd = dragCoeff / dx[0];
        amrex::Vector<amrex::Real> wind{{ux, uy, uz}};
        wind =
            findRefVelocity(verticalSize, x3, device_vel_ht, device_vel_vals);
        // Terrain Drag
        amrex::Real kappa = 0.41;
        amrex::Real ustar =
            std::sqrt(ux * ux + uy * uy) * kappa / std::log((x3 + 0.1) / 0.1);
        amrex::Real Dxz =
            -ustar * ustar * ux / (1e-5 + std::sqrt(ux * ux + uy * uy)) / dx[2];
        amrex::Real Dyz =
            -ustar * ustar * uy / (1e-5 + std::sqrt(ux * ux + uy * uy)) / dx[2];
        // Adjusting Cd for momentum
        amrex::Real CdM = std::min(Cd * 5.0 / (m + 1e-5), 100.0);
        src_term(i, j, k, 0) -=
            (CdM * m * ux * blank(i, j, k) + Dxz * drag(i, j, k) +
             (xdamping + ydamping) * (ux - spongeDensity * wind[0]));
        src_term(i, j, k, 1) -=
            (CdM * m * uy * blank(i, j, k) + Dyz * drag(i, j, k) +
             (xdamping + ydamping) * (uy - spongeDensity * wind[1]));
        src_term(i, j, k, 2) -=
            (CdM * m * uz * blank(i, j, k) +
             (xdamping + ydamping) * (uz - spongeDensity * wind[2]));
    });
}

} // namespace amr_wind::pde::icns