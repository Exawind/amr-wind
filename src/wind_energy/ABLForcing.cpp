#include "ABLForcing.H"

#include "AMReX_ParmParse.H"
#include "AMReX_Gpu.H"

namespace amr_wind {

ABLForcingOld::ABLForcingOld(const SimTime& time)
    : m_time(time)
{
    amrex::ParmParse pp("abl");
    // TODO: Allow forcing at multiple heights
    pp.get("abl_forcing_height", m_forcing_height);

    {
        amrex::ParmParse pp("incflo");
        pp.get("ic_u", m_target_vel[0]);
        pp.get("ic_v", m_target_vel[1]);
        pp.get("ic_w", m_target_vel[2]);
    }

    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
        m_mean_vel[i] = m_target_vel[i];
    }
}

void ABLForcingOld::operator()(
    const amrex::Box& bx,
    const amrex::Array4<amrex::Real>& vel_forces) const
{
    const auto& dt = m_time.deltaT();

    const amrex::Real dudt = (m_target_vel[0] - m_mean_vel[0]) / dt;
    const amrex::Real dvdt = (m_target_vel[1] - m_mean_vel[1]) / dt;

    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        vel_forces(i, j, k, 0) += dudt;
        vel_forces(i, j, k, 1) += dvdt;

        // No forcing in z-direction
    });
}
}
