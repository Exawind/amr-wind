#include "amr-wind/physics/multiphase/DamBreak.H"
#include "amr-wind/physics/multiphase/DamBreakFieldInit.H"
#include "amr-wind/CFDSim.H"
#include "AMReX_ParmParse.H"
#include "amr-wind/fvm/gradient.H"
#include "amr-wind/core/field_ops.H"

namespace amr_wind {

DamBreak::DamBreak(CFDSim& sim)
    : m_velocity(sim.repo().get_field("velocity"))
    , m_levelset(sim.repo().get_field("levelset"))
{
    // Instantiate the DamBreak field initializer
    m_field_init.reset(new DamBreakFieldInit());
}

/** Initialize the velocity and levelset fields at the beginning of the
 *  simulation.
 *
 *  \sa amr_wind::DamBreakFieldInit
 */
void DamBreak::initialize_fields(int level, const amrex::Geometry& geom)
{
    auto& velocity = m_velocity(level);
    auto& levelset = m_levelset(level);

    for (amrex::MFIter mfi(levelset); mfi.isValid(); ++mfi) {
        const auto& vbx = mfi.validbox();

        (*m_field_init)(vbx, geom, velocity.array(mfi), levelset.array(mfi));
    }
}

} // namespace amr_wind