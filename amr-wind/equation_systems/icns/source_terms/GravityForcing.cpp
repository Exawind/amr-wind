#include "amr-wind/equation_systems/icns/source_terms/GravityForcing.H"
#include "amr-wind/CFDSim.H"
#include "amr-wind/core/FieldUtils.H"

#include "AMReX_ParmParse.H"

namespace amr_wind::pde::icns {

/** Gravity Forcing source term
 *
 *  Reads in the following parameters from `incflo` namespace:
 *
 *  - `gravity` acceleration due to gravity (m/s)
 */
GravityForcing::GravityForcing(const CFDSim& sim)
{
    amrex::ParmParse pp("incflo");
    pp.queryarr("gravity", m_gravity);
    pp.query("density", m_rho0_const);

    // Get density fields
    m_rho = &(sim.repo().get_field("density"));

    // Check if perturbational pressure desired
    amrex::ParmParse pp_icns("ICNS");
    pp_icns.query("use_perturb_pressure", is_pptb);
    // Determine if rho0 field exists
    is_rho0 = sim.repo().field_exists("reference_density");
    if (is_rho0) {
        m_rho0 = &(sim.repo().get_field("reference_density"));
    } else {
        // Point to existing field, won't be used
        m_rho0 = m_rho;
    }
}

GravityForcing::~GravityForcing() = default;

/** Add the Gravity source term to the forcing array
 *
 *  @param lev AMR level
 *  @param mfi multiFab index
 *  @param bx Box to operate on
 *  @param FieldState field
 *  @param vel_forces Forcing source term
 */
void GravityForcing::operator()(
    const int lev,
    const amrex::MFIter& mfi,
    const amrex::Box& bx,
    const FieldState fstate,
    const amrex::Array4<amrex::Real>& vel_forces) const
{
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> gravity{
        {m_gravity[0], m_gravity[1], m_gravity[2]}};

    const auto& rho_arr =
        ((*m_rho).state(field_impl::phi_state(fstate)))(lev).const_array(mfi);
    const auto& rho0_arr = (*m_rho0)(lev).const_array(mfi);
    const bool ir0 = is_rho0;
    const bool ipt = is_pptb;
    const amrex::Real mr0c = m_rho0_const;

    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        const amrex::Real factor =
            (!ipt ? 1.0
                  : 1.0 - (ir0 ? rho0_arr(i, j, k) : mr0c) / rho_arr(i, j, k));

        vel_forces(i, j, k, 0) += gravity[0] * factor;
        vel_forces(i, j, k, 1) += gravity[1] * factor;
        vel_forces(i, j, k, 2) += gravity[2] * factor;
    });
}

} // namespace amr_wind::pde::icns
