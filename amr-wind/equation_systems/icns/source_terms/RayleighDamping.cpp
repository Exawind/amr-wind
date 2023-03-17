#include "amr-wind/equation_systems/icns/source_terms/RayleighDamping.H"
#include "amr-wind/CFDSim.H"
#include "amr-wind/utilities/trig_ops.H"

#include "AMReX_ParmParse.H"
#include "AMReX_Gpu.H"

namespace amr_wind::pde::icns {

/** Body Force
 */
RayleighDamping::RayleighDamping(const CFDSim& sim)
    : m_mesh(sim.mesh()), m_velocity(sim.repo().get_field("velocity"))
{
    // Rea the Rayleigh Damping Layer parameters
    amrex::ParmParse pp("RayleighDamping");
    pp.query("time_scale", m_tau);
    pp.query("damping_length", m_dRD);
    pp.query("full_damping_length", m_dFull);
    pp.getarr("reference_velocity", m_ref_vel);
}

RayleighDamping::~RayleighDamping() = default;

void RayleighDamping::operator()(
    const int lev,
    const amrex::MFIter& mfi,
    const amrex::Box& bx,
    const FieldState fstate,
    const amrex::Array4<amrex::Real>& src_term) const
{

    const auto& problo = m_mesh.Geom(lev).ProbLoArray();
    const auto& probhi = m_mesh.Geom(lev).ProbHiArray();
    const auto& dx = m_mesh.Geom(lev).CellSizeArray();

    const amrex::Real tau = m_tau;
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> ref_vel{
        {m_ref_vel[0], m_ref_vel[1], m_ref_vel[2]}};
    const auto& vel =
        m_velocity.state(field_impl::dof_state(fstate))(lev).const_array(mfi);

    // Constants used to determine the fringe region coefficient
    const amrex::Real dRD = m_dRD;
    const amrex::Real dFull = m_dFull;

    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        amrex::Real coeff = 0.0; // Need to work on defining this
        const amrex::Real z = problo[2] + (k + 0.5) * dx[2];

        if (probhi[2] - z > dRD) {
            coeff = 0.0;
        } else if (probhi[2] - z > dFull) {
            coeff =
                0.5 * std::cos((probhi[2] - dFull - z) / (dRD - dFull)) + 0.5;
        } else {
            coeff = 1.0;
        }
        src_term(i, j, k, 0) += coeff * (vel(i, j, k, 0) - ref_vel[0]) / tau;
        src_term(i, j, k, 1) += coeff * (vel(i, j, k, 1) - ref_vel[1]) / tau;
        src_term(i, j, k, 2) += coeff * (vel(i, j, k, 2) - ref_vel[2]) / tau;
    });
}

} // namespace amr_wind::pde::icns