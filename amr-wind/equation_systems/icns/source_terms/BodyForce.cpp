#include "amr-wind/equation_systems/icns/source_terms/BodyForce.H"
#include "amr-wind/CFDSim.H"
#include "amr-wind/utilities/trig_ops.H"

#include "AMReX_ParmParse.H"
#include "AMReX_Gpu.H"

namespace amr_wind::pde::icns {

/** Body Force
 */
BodyForce::BodyForce(const CFDSim& sim) : m_time(sim.time())
{

    // Read the geostrophic wind speed vector (in m/s)
    amrex::ParmParse pp("BodyForce");
    pp.query("type", m_type);
    m_type = amrex::toLower(m_type);
    pp.getarr("magnitude", m_body_force);
    if (m_type == "oscillatory") {
        pp.get("angular_frequency", m_omega);
    }
}

BodyForce::~BodyForce() = default;

void BodyForce::operator()(
    const int /*lev*/,
    const amrex::MFIter& /*mfi*/,
    const amrex::Box& bx,
    const FieldState /*fstate*/,
    const amrex::Array4<amrex::Real>& src_term) const
{
    const auto& time = m_time.current_time();
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> forcing{
        {m_body_force[0], m_body_force[1], m_body_force[2]}};

    amrex::Real coeff =
        (m_type == "oscillatory") ? std::cos(m_omega * time) : 1.0;

    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        src_term(i, j, k, 0) += coeff * forcing[0];
        src_term(i, j, k, 1) += coeff * forcing[1];
        src_term(i, j, k, 2) += coeff * forcing[2];
    });
}

} // namespace amr_wind::pde::icns
