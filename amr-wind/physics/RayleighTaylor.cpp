#include "amr-wind/physics/RayleighTaylor.H"
#include "amr-wind/physics/RayleighTaylorFieldInit.H"
#include "amr-wind/CFDSim.H"

namespace amr_wind {

RayleighTaylor::RayleighTaylor(const CFDSim& sim)
    : m_velocity(sim.repo().get_field("velocity"))
    , m_density(sim.repo().get_field("density"))
{
    // Instantiate the RayleighTaylor field initializer
    m_field_init.reset(new RayleighTaylorFieldInit());
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

    velocity.setVal(0.0);

    for (amrex::MFIter mfi(density); mfi.isValid(); ++mfi) {
        const auto& vbx = mfi.validbox();

        (*m_field_init)(
            vbx, geom, density.array(mfi));
    }
}

} // namespace amr_wind
