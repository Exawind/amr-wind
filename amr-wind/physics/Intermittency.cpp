#include "amr-wind/physics/Intermittency.H"
#include "amr-wind/CFDSim.H"
#include "AMReX_ParmParse.H"
#include "amr-wind/utilities/trig_ops.H"

namespace amr_wind {

Intermittency::Intermittency(const CFDSim& sim)
    : m_intermittency(sim.repo().declare_field("intermittency",1,0,1))
{
    amrex::ParmParse pp("incflo");
    pp.query("gamma_intermittency", m_gamma);
}

/** Initialize the intermittency field at the beginning of the
 *  simulation.
 */
void Intermittency::initialize_fields(
    int level, const amrex::Geometry& geom)
{
    using namespace utils;

    auto& gamma = m_intermittency(level);

    gamma.setVal(m_gamma);
}

} // namespace amr_wind
