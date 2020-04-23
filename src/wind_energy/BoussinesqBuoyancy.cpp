#include "BoussinesqBuoyancy.H"
#include "CFDSim.H"
#include "FieldUtils.H"

#include "AMReX_ParmParse.H"

namespace amr_wind {

/** Boussinesq buoyancy source term for ABL simulations
 *
 *  Reads in the following parameters from `abl` namespace:
 *
 *  - `reference_temperature` (Mandatory) temperature (`T0`) in Kelvin
 *  - `thermal_expansion_coeff` Optional, default = `1.0 / T0`
 *  - `gravity` acceleration due to gravity (m/s)
 */
BoussinesqBuoyancyOld::BoussinesqBuoyancyOld()
{
    // fixme: do we want to use abl namespace if it can be used by other physics?
    amrex::ParmParse pp("abl");

    pp.get("reference_temperature", m_ref_theta);

    if (pp.contains("thermal_expansion_coeff")) {
        pp.get("thermal_expansion_coeff", m_beta);
    } else {
        m_beta = 1.0 / m_ref_theta;
    }

    // FIXME: gravity in `incflo` namespace
    {
        amrex::ParmParse pp("incflo");
        pp.queryarr("gravity", m_gravity);
    }
}

/** Add the Boussinesq source term to the forcing array
 *
 *  @param bx Box to operate on
 *  @param scalars Transported scalars (first one is assumed to be temperature)
 *  @param vel_forces Forcing source term where Boussinesq term is added
 */
void BoussinesqBuoyancyOld::operator()(
    const amrex::Box& bx,
    const amrex::Array4<const amrex::Real>& scalars,
    const amrex::Array4<amrex::Real>& vel_forces) const
{
    const amrex::Real T0 = m_ref_theta;
    const amrex::Real beta = m_beta;
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> gravity{
        m_gravity[0], m_gravity[1], m_gravity[2]};

    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        const amrex::Real T = scalars(i, j, k, 0);
        const amrex::Real fac = beta * (T0 - T);

        vel_forces(i, j, k, 0) += gravity[0] * fac;
        vel_forces(i, j, k, 1) += gravity[1] * fac;
        vel_forces(i, j, k, 2) += gravity[2] * fac;
    });
}

namespace pde {
namespace icns {

BoussinesqBuoyancy::BoussinesqBuoyancy(const CFDSim& sim)
    : m_temperature(sim.repo().get_field("temperature"))
{
    amrex::ParmParse pp(this->identifier());

    pp.get("reference_temperature", m_ref_theta);

    if (pp.contains("thermal_expansion_coeff")) {
        pp.get("thermal_expansion_coeff", m_beta);
    } else {
        m_beta = 1.0 / m_ref_theta;
    }

    // FIXME: gravity in `incflo` namespace
    {
        amrex::ParmParse pp("incflo");
        pp.queryarr("gravity", m_gravity);
    }
}

BoussinesqBuoyancy::~BoussinesqBuoyancy() = default;

void BoussinesqBuoyancy::operator()(
    const int lev,
    const amrex::MFIter& mfi,
    const amrex::Box& bx,
    const FieldState fstate,
    const amrex::Array4<amrex::Real>& src_term) const
{
    const amrex::Real T0 = m_ref_theta;
    const amrex::Real beta = m_beta;
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> gravity{
        m_gravity[0], m_gravity[1], m_gravity[2]};

    const auto& temp =
        m_temperature.state(field_impl::phi_state(fstate))(lev).const_array(
            mfi);
    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        const amrex::Real T = temp(i, j, k, 0);
        const amrex::Real fac = beta * (T0 - T);

        src_term(i, j, k, 0) += gravity[0] * fac;
        src_term(i, j, k, 1) += gravity[1] * fac;
        src_term(i, j, k, 2) += gravity[2] * fac;
    });
}

} // namespace icns
} // namespace pde
} // namespace amr_wind
