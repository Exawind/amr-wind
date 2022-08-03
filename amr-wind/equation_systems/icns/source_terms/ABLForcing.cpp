#include "amr-wind/equation_systems/icns/source_terms/ABLForcing.H"
#include "amr-wind/CFDSim.H"
#include "amr-wind/wind_energy/ABL.H"
#include "amr-wind/utilities/trig_ops.H"

#include "AMReX_ParmParse.H"
#include "AMReX_Gpu.H"

namespace amr_wind {
namespace pde {
namespace icns {

ABLForcing::ABLForcing(const CFDSim& sim) : m_time(sim.time())
{
    const auto& abl = sim.physics_manager().get<amr_wind::ABL>();
    abl.register_forcing_term(this);
    abl.abl_statistics().register_forcing_term(this);

    amrex::ParmParse pp_abl(identifier());
    // TODO: Allow forcing at multiple heights
    pp_abl.get("abl_forcing_height", m_forcing_height);
    amrex::ParmParse pp_incflo("incflo");

    if (pp_abl.contains("velocity_timetable")) {
        // MICHAEL: interpolate target velocity
        m_speed = m_vel_speed[0];
        m_dir = m_vel_dir[0];
        for (int it = 0; it < num_vel_pts - 1; ++it) {
            if ((m_time > m_vel_time[it]) && (m_time <= m_vel_time[it + 1])) {
                const amrex::Real slope_spd =
                    (m_vel_speed[it + 1] - m_vel_speed[it]) / (m_vel_time[it + 1] - m_vel_time[it]);
                const amrex::Real slope_dir =
                    (m_vel_dir[it + 1] - m_vel_dir[it]) / (m_vel_time[it + 1] - m_vel_time[it]);
                m_speed = m_vel_speed[it] + (m_time - m_vel_time[it]) * slope_spd;
                m_dir = m_vel_dir[it] + (m_time - m_vel_time[it]) * slope_dir;
            }
        }

        // Convert to Cartesian coordinates
        m_dir_rad = utils::radians(m_dir);
        m_target_vel[0] = m_speed * std::cos(m_dir_rad);
        m_target_vel[1] = m_speed * std::sin(m_dir_rad);

    } else {
        pp_incflo.getarr("velocity", m_target_vel);
    }

    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
        m_mean_vel[i] = m_target_vel[i];
    }
}

ABLForcing::~ABLForcing() = default;

void ABLForcing::operator()(
    const int /*lev*/,
    const amrex::MFIter& /*mfi*/,
    const amrex::Box& bx,
    const FieldState /*fstate*/,
    const amrex::Array4<amrex::Real>& src_term) const
{

    const amrex::Real dudt = m_abl_forcing[0];
    const amrex::Real dvdt = m_abl_forcing[1];

    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        src_term(i, j, k, 0) += dudt;
        src_term(i, j, k, 1) += dvdt;

        // No forcing in z-direction
    });
}

} // namespace icns
} // namespace pde
} // namespace amr_wind
