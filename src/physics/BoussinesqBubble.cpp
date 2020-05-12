#include "BoussinesqBubble.H"
#include "BoussinesqBubbleFieldInit.H"
#include "CFDSim.H"
#include "ABLWallFunction.H"

namespace amr_wind {

BoussinesqBubble::BoussinesqBubble(CFDSim& sim)
    : m_velocity(sim.repo().get_field("velocity"))
    , m_density(sim.repo().get_field("density"))
{
    // Register temperature equation
    auto& teqn = sim.pde_manager().register_transport_pde("Temperature");

    // Defer getting temperature field until PDE has been registered
    m_temperature = &(teqn.fields().field);

    // Instantiate the BoussinesqBubble field initializer
    m_field_init.reset(new BoussinesqBubbleFieldInit());
}

/** Initialize the velocity and temperature fields at the beginning of the
 *  simulation.
 *
 *  \sa amr_wind::BoussinesqBubbleFieldInit
 */
void BoussinesqBubble::initialize_fields(
    int level,
    const amrex::Geometry& geom)
{
    auto& velocity = m_velocity(level);
    auto& density = m_density(level);
    auto& scalars = (*m_temperature)(level);

    for (amrex::MFIter mfi(density); mfi.isValid(); ++mfi) {
        const auto& vbx = mfi.validbox();

        (*m_field_init)(
            vbx, geom, velocity.array(mfi), density.array(mfi),
            scalars.array(mfi));
    }
}

} // namespace amr_wind
