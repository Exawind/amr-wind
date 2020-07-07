#include "amr-wind/equation_systems/icns/source_terms/BodyForce.H"
#include "amr-wind/CFDSim.H"
#include "amr-wind/utilities/trig_ops.H"

#include "AMReX_ParmParse.H"
#include "AMReX_Gpu.H"

namespace amr_wind {
namespace pde {
namespace icns {

/** Body Force
 */
BodyForce::BodyForce(const CFDSim& sim)
    : m_time(sim.time()), m_density(sim.repo().get_field("density"))
{

    {
        // Read the geostrophic wind speed vector (in m/s)
        amrex::ParmParse pp("BodyForce");
        pp.query("type", m_type);
        pp.queryarr("magnitude", m_body_force);
        if (m_type == "Oscillatory" || m_type == "oscillatory") {
            pp.query("AngularFrequency", m_omega);
        }
    }
}

BodyForce::~BodyForce() = default;

void BodyForce::operator()(
    const int lev,
    const amrex::MFIter& mfi,
    const amrex::Box& bx,
    const FieldState fstate,
    const amrex::Array4<amrex::Real>& src_term) const
{

    const auto& time = m_time.current_time();
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> forcing{
        {m_body_force[0], m_body_force[1], m_body_force[2]}};

    amrex::Real coeff;

    if (m_type == "constant" || m_type == "Constant") {
        coeff = 1.;
    } else if (m_type == "Oscillatory" || m_type == "oscillatory") {
        coeff = std::sin(m_omega * time);
    }
    const auto& rho = m_density.state(fstate)(lev).const_array(mfi);
    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        src_term(i, j, k, 0) += rho(i, j, k) * coeff * forcing[0];
        src_term(i, j, k, 1) += rho(i, j, k) * coeff * forcing[1];
        src_term(i, j, k, 2) += rho(i, j, k) * coeff * forcing[2];
    });
}

} // namespace icns
} // namespace pde
} // namespace amr_wind
