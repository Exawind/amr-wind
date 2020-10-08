#include "amr-wind/physics/VortexPatch.H"
#include "amr-wind/physics/VortexPatchFieldInit.H"
#include "amr-wind/CFDSim.H"

namespace amr_wind {

VortexPatch::VortexPatch(CFDSim& sim)
    : m_velocity(sim.repo().get_field("velocity"))
    , m_density(sim.repo().get_field("density"))
{
    // Register levelset equation
    auto& levelset_eqn = sim.pde_manager().register_transport_pde("Levelset");

    // Defer getting levelset field until PDE has been registered
    m_levelset = &(levelset_eqn.fields().field);

    // Instantiate the VortexPatch field initializer
    m_field_init.reset(new VortexPatchFieldInit());
}

/** Initialize the velocity and levelset fields at the beginning of the
 *  simulation.
 *
 *  \sa amr_wind::VortexPatchFieldInit
 */
void VortexPatch::initialize_fields(
    int level,
    const amrex::Geometry& geom)
{
    auto& velocity = m_velocity(level);
    auto& density = m_density(level);
    auto& levelset = (*m_levelset)(level);

    for (amrex::MFIter mfi(density); mfi.isValid(); ++mfi) {
        const auto& vbx = mfi.validbox();

        (*m_field_init)(
            vbx, geom, velocity.array(mfi), density.array(mfi),
            levelset.array(mfi));
    }
}

} // namespace amr_wind