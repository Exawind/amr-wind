#include <memory>

#include "amr-wind/physics/RayleighTaylor.H"
#include "amr-wind/physics/RayleighTaylorFieldInit.H"
#include "amr-wind/CFDSim.H"
#include "AMReX_ParmParse.H"

namespace amr_wind {

RayleighTaylor::RayleighTaylor(const CFDSim& sim)
    : m_velocity(sim.repo().get_field("velocity"))
    , m_density(sim.repo().get_field("density"))
    , m_rho0(sim.repo().declare_field("rho0", 1, 0, 1))
    , m_field_init(std::make_unique<RayleighTaylorFieldInit>())
{
    // Check for non-default reference density
    amrex::ParmParse pp("incflo");
    pp.query("density", m_rho0_const);
}

/** Initialize the velocity and density fields at the beginning of the
 *  simulation.
 *
 *  \sa amr_wind::RayleighTaylorFieldInit
 */
void RayleighTaylor::initialize_fields(int level, const amrex::Geometry& geom)
{
    auto& velocity = m_velocity(level);
    auto& density = m_density(level);
    auto& rho0 = m_rho0(level);

    velocity.setVal(0.0);

    // Initialize reference density field
    rho0.setVal(m_rho0_const);

    for (amrex::MFIter mfi(density); mfi.isValid(); ++mfi) {
        const auto& vbx = mfi.validbox();

        (*m_field_init)(vbx, geom, density.array(mfi));
    }
}

void RayleighTaylor::post_regrid_actions() { m_rho0.setVal(m_rho0_const); }

} // namespace amr_wind
