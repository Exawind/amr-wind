#include "FreeStream.H"
#include "CFDSim.H"
#include "AMReX_ParmParse.H"
#include "trig_ops.H"

namespace amr_wind {

FreeStream::FreeStream(const CFDSim& sim)
    : m_velocity(sim.repo().get_field("velocity"))
    , m_density(sim.repo().get_field("density"))
{
    amrex::ParmParse pp("FreeStream");

    {
        amrex::ParmParse pp("incflo");
        pp.query("ro_0", m_rho);
        pp.query("ic_u", m_u);
        pp.query("ic_v", m_v);
        pp.query("ic_w", m_w);
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
    velocity.setVal(m_u, 0, 1);
    velocity.setVal(m_v, 1, 1);
    velocity.setVal(m_w, 2, 1);

}

} // namespace amr_wind
