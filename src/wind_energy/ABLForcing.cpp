#include "ABLForcing.H"

#include "AMReX_ParmParse.H"

namespace amr_wind {

ABLForcing::ABLForcing(const SimTime& time)
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
}

void ABLForcing::operator()(
    const amrex::Box& vbx,
    const amrex::Array4<amrex::Real>& velocity,
    amrex::Array4<amrex::Real>& vel_forces) const
{
    const auto& dt = m_time.deltaT();
}

}
