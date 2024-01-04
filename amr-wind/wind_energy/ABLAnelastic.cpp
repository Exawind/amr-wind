#include "amr-wind/wind_energy/ABL.H"

#include "AMReX_ParmParse.H"
#include "AMReX_Print.H"

namespace amr_wind {

ABLAnelastic::ABLAnelastic(CFDSim& sim) : m_mesh(sim.mesh())
{
    amrex::ParmParse pp("ABL");
    pp.query("anelastic", m_is_anelastic);
}

void ABLAnelastic::post_init_actions()
{
    if (!m_is_anelastic) {
        return;
    }
    initialize_data();
}
void ABLAnelastic::initialize_data()
{
    amrex::Print() << "hello!!" << std::endl;
}
} // namespace amr_wind
