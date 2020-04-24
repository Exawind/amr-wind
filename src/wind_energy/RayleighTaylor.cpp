#include "RayleighTaylor.H"
#include "RayleighTaylorFieldInit.H"
#include "CFDSim.H"

namespace amr_wind {

RayleighTaylor::RayleighTaylor(const CFDSim& sim)
    : m_velocity(sim.repo().get_field("velocity"))
    , m_density(sim.repo().get_field("density"))
    , m_temperature(sim.repo().get_field("temperature"))
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
    auto& scalars = m_temperature(level);

    for (amrex::MFIter mfi(density); mfi.isValid(); ++mfi) {
        const auto& vbx = mfi.validbox();

        (*m_field_init)(
            vbx, geom, velocity.array(mfi), density.array(mfi),
            scalars.array(mfi));
    }
}

} // namespace amr_wind
