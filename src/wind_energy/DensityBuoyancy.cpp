//
//  DensityBuoyancy.cpp
//  amr-wind
//

#include "DensityBuoyancy.H"

#include "AMReX_ParmParse.H"

namespace amr_wind {

/** Density buoyancy source term
 *
 *  Reads in the following parameters from `incflo` namespace:
 *
 *  - `gravity` acceleration due to gravity (m/s)
 *  - `ro_0` reference density (kg/m^3)
 */
DensityBuoyancy::DensityBuoyancy(const SimTime&, const FieldRepo& repo_in)
    : m_density(repo_in.get_field("density"))
{
    amrex::ParmParse pp("incflo");
    pp.queryarr("gravity", m_gravity);
    pp.query("ro_0", rho_0);
}

/** Add the Density source term to the forcing array
 *
 *  @param bx Box to operate on
 *  @param density scalar
 *  @param vel_forces Forcing source term where buoyancy from density is added
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

} // namespace amr_wind
