#include "amr-wind/physics/FreeStream.H"
#include "amr-wind/CFDSim.H"
#include "AMReX_ParmParse.H"
#include "amr-wind/utilities/trig_ops.H"

namespace amr_wind {

FreeStream::FreeStream(const CFDSim& sim)
    : m_velocity(sim.repo().get_field("velocity"))
    , m_density(sim.repo().get_field("density"))
{

    // no freestream inputs yet
    // amrex::ParmParse pp(this->identifier());

    // fixme rho,u,v,w in incflo class
    {
        amrex::ParmParse pp("incflo");
        pp.query("density", m_rho);
        pp.queryarr("velocity",m_vel);

    }
}

/** Initialize the velocity and density fields at the beginning of the
 *  simulation.
 */
void FreeStream::initialize_fields(
    int level,
    const amrex::Geometry& /* geom */)
{
    auto& velocity = m_velocity(level);
    auto& density = m_density(level);

    density.setVal(m_rho);
    velocity.setVal(m_vel[0], 0, 1);
    velocity.setVal(m_vel[1], 1, 1);
    velocity.setVal(m_vel[2], 2, 1);

}

} // namespace amr_wind
