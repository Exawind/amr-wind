#include "amr-wind/equation_systems/icns/source_terms/GravityForcing.H"
#include "amr-wind/CFDSim.H"
#include "amr-wind/core/FieldUtils.H"

#include "AMReX_ParmParse.H"

namespace amr_wind {
namespace pde {
namespace icns {

/** Gracity Forcing source term
 *
 *  Reads in the following parameters from `incflo` namespace:
 *
 *  - `gravity` acceleration due to gravity (m/s)
 */
GravityForcing::GravityForcing(const CFDSim& sim)
    : m_density(sim.repo().get_field("density"))
{
    {
        amrex::ParmParse pp("incflo");
        pp.queryarr("gravity", m_gravity);
    }
}

GravityForcing::~GravityForcing() = default;

/** Add the Gravity source term to the forcing array
 *
 *  @param lev AMR level
 *  @param mfi multiFab index
 *  @param bx Box to operate on
 *  @param FieldState field
 *  @param vel_forces Forcing source term, activated when density varies from
 * rho_0
 */
void GravityForcing::operator()(
    const int,
    const amrex::MFIter&,
    const amrex::Box& bx,
    const FieldState,
    const amrex::Array4<amrex::Real>& vel_forces) const
{
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> gravity{
        {m_gravity[0], m_gravity[1], m_gravity[2]}};

    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        vel_forces(i, j, k, 0) += gravity[0];
        vel_forces(i, j, k, 1) += gravity[1];
        vel_forces(i, j, k, 2) += gravity[2];
    });
}

} // namespace icns
} // namespace pde
} // namespace amr_wind