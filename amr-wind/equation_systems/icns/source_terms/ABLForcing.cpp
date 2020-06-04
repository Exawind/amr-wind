#include "amr-wind/equation_systems/icns/source_terms/ABLForcing.H"
#include "amr-wind/utilities/PlaneAveraging.H"
#include "amr-wind/CFDSim.H"
#include "amr-wind/wind_energy/ABL.H"

#include "AMReX_ParmParse.H"
#include "AMReX_Gpu.H"

namespace amr_wind {
namespace pde {
namespace icns {

ABLForcing::ABLForcing(const CFDSim& sim) : m_time(sim.time())
{
    const auto& abl = sim.physics_manager().get<amr_wind::ABL>();
    abl.register_forcing_term(this);

    amrex::ParmParse pp_abl(identifier());
    // TODO: Allow forcing at multiple heights
    pp_abl.get("abl_forcing_height", m_forcing_height);
    amrex::ParmParse pp_incflo("incflo");
    pp_incflo.getarr("velocity", m_target_vel);

    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
        m_mean_vel[i] = m_target_vel[i];
    }
}

ABLForcing::~ABLForcing() = default;

void ABLForcing::operator()(
    const int,
    const amrex::MFIter&,
    const amrex::Box& bx,
    const FieldState,
    const amrex::Array4<amrex::Real>& src_term) const
{
    const auto& dt = m_time.deltaT();

    const amrex::Real dudt = (m_target_vel[0] - m_mean_vel[0]) / dt;
    const amrex::Real dvdt = (m_target_vel[1] - m_mean_vel[1]) / dt;

    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        src_term(i, j, k, 0) += dudt;
        src_term(i, j, k, 1) += dvdt;

        // No forcing in z-direction
    });
}

} // namespace icns
} // namespace pde
} // namespace amr_wind
