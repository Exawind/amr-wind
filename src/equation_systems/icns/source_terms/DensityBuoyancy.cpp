//
//  DensityBuoyancy.cpp
//  amr-wind
//

#include "DensityBuoyancy.H"
#include "CFDSim.H"
#include "FieldUtils.H"

#include "AMReX_ParmParse.H"

namespace amr_wind {
namespace pde {
namespace icns {

/** Density based buoyancy source term
 *
 *  Reads in the following parameters from `incflo` namespace:
 *
 *  - `gravity` acceleration due to gravity (m/s)
 *  - `reference density` Optional, default = `1.0`
 */
DensityBuoyancy::DensityBuoyancy(const CFDSim& sim)
    : m_density(sim.repo().get_field("density"))
{
    amrex::ParmParse pp(this->identifier());

    // FIXME: gravity in `incflo` namespace
    {
        amrex::ParmParse pp("incflo");
        pp.queryarr("gravity", m_gravity);
        pp.query("ro_0", rho_0);
    }
}

DensityBuoyancy::~DensityBuoyancy() = default;

/** Add the Boussinesq source term to the forcing array
 *
 *  @param lev AMR level
 *  @param mfi multiFab index
 *  @param bx Box to operate on
 *  @param fstate field state
 *  @param vel_forces Forcing source term, activated when density varies from
 * rho_0
 */
void DensityBuoyancy::operator()(
    const int lev,
    const amrex::MFIter& mfi,
    const amrex::Box& bx,
    const FieldState fstate,
    const amrex::Array4<amrex::Real>& vel_forces) const
{
    const amrex::Real density_0 = rho_0;
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> gravity{
        m_gravity[0], m_gravity[1], m_gravity[2]};

    FieldState den_state = field_impl::phi_state(fstate);
    const auto& density = m_density.state(den_state)(lev).const_array(mfi);

    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {

        const amrex::Real fac = 1.0 - density_0 / density(i, j, k);

        vel_forces(i, j, k, 0) += gravity[0] * fac;
        vel_forces(i, j, k, 1) += gravity[1] * fac;
        vel_forces(i, j, k, 2) += gravity[2] * fac;
    });
}

} // namespace icns
} // namespace pde
} // namespace amr_wind
