#include "amr-wind/equation_systems/icns/source_terms/ABLMeanBoussinesq.H"
#include "amr-wind/CFDSim.H"
#include "amr-wind/core/FieldUtils.H"

#include "AMReX_ParmParse.H"

namespace amr_wind {
namespace pde {
namespace icns {

/** Boussinesq buoyancy source term for ABL simulations
 *
 *  Reads in the following parameters from `ABLMeanBoussinesq` namespace:
 *
 *  - `reference_temperature` (Mandatory) temperature (`T0`) in Kelvin
 *  - `thermal_expansion_coeff` Optional, default = `1.0 / T0`
 *  - `gravity` acceleration due to gravity (m/s)
 */
ABLMeanBoussinesq::ABLMeanBoussinesq(const CFDSim& sim)
    : m_mesh(sim.mesh())
{
    amrex::ParmParse pp_boussinesq_buoyancy("BoussinesqBuoyancy");
    pp_boussinesq_buoyancy.get("reference_temperature", m_ref_theta);

    if (pp_boussinesq_buoyancy.contains("thermal_expansion_coeff")) {
        pp_boussinesq_buoyancy.get("thermal_expansion_coeff", m_beta);
    } else {
        m_beta = 1.0 / m_ref_theta;
    }

    // FIXME: gravity in `incflo` namespace
    amrex::ParmParse pp_incflo("incflo");
    pp_incflo.queryarr("gravity", m_gravity);
}

ABLMeanBoussinesq::~ABLMeanBoussinesq() = default;

void ABLMeanBoussinesq::operator()(
    const int lev,
    const amrex::MFIter&,
    const amrex::Box& bx,
    const FieldState,
    const amrex::Array4<amrex::Real>& src_term) const
{
    const auto& problo = m_mesh.Geom(lev).ProbLoArray();
    const auto& dx = m_mesh.Geom(lev).CellSizeArray();
    const amrex::Real T0 = m_ref_theta;
    const amrex::Real beta = m_beta;
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> gravity{{
        m_gravity[0], m_gravity[1], m_gravity[2]}};

    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        amrex::Real temp = 308.75;
        const amrex::Real z = problo[2] + (k + 0.5) * dx[2];
        if (z > 750.0) {
            temp = 308.0 + 0.003 * (z - 750.0);
        } else if (z > 650.0) {
            temp = 300.0 + 0.08 * (z - 650.0);
        } else {
            temp = 300.0;
        }

        const amrex::Real fac = beta * (temp - T0);
        src_term(i, j, k, 0) += gravity[0] * fac;
        src_term(i, j, k, 1) += gravity[1] * fac;
        src_term(i, j, k, 2) += gravity[2] * fac;
    });
}

} // namespace icns
} // namespace pde
} // namespace amr_wind
