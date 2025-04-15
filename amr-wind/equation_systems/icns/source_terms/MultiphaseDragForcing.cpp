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
    , m_vof(sim.repo().get_field("vof"))
    , m_density(sim.repo().get_field("density"))
{
    amrex::ParmParse pp(identifier());
    pp.query("wave_roughness", m_wave_roughness);
    pp.query("density_ratio_limit", m_density_ratio_limit);
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
    const auto& vof = m_vof(lev).const_array(mfi);
    const auto& rho = m_density(lev).const_array(mfi);
    const auto& geom = m_mesh.Geom(lev);
    const auto& dx = geom.CellSizeArray();
    const auto& dt = m_time.delta_t();
    const amrex::Real kappa = 0.41;
    const amrex::Real z0 = m_wave_roughness;
    const amrex::Real f_rho = m_density_ratio_limit;
    const auto tiny = std::numeric_limits<amrex::Real>::epsilon();
    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        const amrex::IntVect iv{i, j, k};
        amrex::IntVect iv_gas{iv};
        amrex::IntVect iv_liquid{iv};
        int add_drag{0};
        if (vof(i, j, k) < 0.5 && rho(i, j, k) * f_rho < rho(i, j, k - 1)) {
            // Current cell is mostly gas, density ratio indicates drag should
            // be added
            add_drag = 1;
            iv_liquid[2] -= 1;
        }
        if (vof(i, j, k + 1) < 0.5 && rho(i, j, k + 1) * f_rho < rho(i, j, k)) {
            // Above cell is mostly gas, density ratio indicates drag should be
            // subtracted
            add_drag = -1;
            iv_gas[2] += 1;
        }
        const amrex::IntVect iv_above_gas{iv_gas[0], iv_gas[1], iv_gas[2] + 1};

        // Strain rate uses gas for flow velocity, liquid for wall velocity
        const amrex::Real ux1 = vel(iv_gas, 0);
        const amrex::Real uy1 = vel(iv_gas, 1);
        amrex::Real Dxz = 0.0;
        amrex::Real Dyz = 0.0;
        amrex::Real bc_forcing_x = 0;
        amrex::Real bc_forcing_y = 0;

        const amrex::Real wall_u = vel(iv_liquid, 0);
        const amrex::Real wall_v = vel(iv_liquid, 1);
        // Relative velocities for calculating shear
        const amrex::Real ux1r = ux1 - wall_u;
        const amrex::Real uy1r = uy1 - wall_v;
        const amrex::Real ux2r = vel(iv_above_gas, 0) - wall_u;
        const amrex::Real uy2r = vel(iv_above_gas, 1) - wall_v;
        const amrex::Real ustar = viscous_drag_calculations(
            Dxz, Dyz, ux1r, uy1r, ux2r, uy2r, z0, dx[2], kappa, tiny);
        const amrex::Real uTarget = ustar / kappa * std::log(0.5 * dx[2] / z0);
        const amrex::Real uxTarget =
            uTarget * ux2r / (tiny + std::sqrt(ux2r * ux2r + uy2r * uy2r));
        const amrex::Real uyTarget =
            uTarget * uy2r / (tiny + std::sqrt(ux2r * ux2r + uy2r * uy2r));
        // BC forcing pushes nonrelative velocity toward target velocity
        bc_forcing_x = -(uxTarget - ux1) / dt;
        bc_forcing_y = -(uyTarget - uy1) / dt;

        // Drag force is calculated in gas phase and is equal and opposite in
        // neighboring cell. Because source term is kinematic and local density
        // is multiplied later, need to apply density ratio here
        const amrex::Real r_rho = rho(iv_gas) / rho(iv);

        src_term(i, j, k, 0) -=
            static_cast<amrex::Real>(add_drag) * r_rho * (Dxz + bc_forcing_x);
        src_term(i, j, k, 1) -=
            static_cast<amrex::Real>(add_drag) * r_rho * (Dyz + bc_forcing_y);
    });
}

} // namespace amr_wind::pde::icns
