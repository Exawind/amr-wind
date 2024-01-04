#include "amr-wind/wind_energy/ABL.H"

namespace amr_wind {

ABLAnelastic::ABLAnelastic(CFDSim& sim) : m_mesh(sim.mesh())
{
    amrex::ParmParse pp("ABL");
    pp.query("anelastic", m_is_anelastic);
    pp.query("anelastic_axis", m_axis);
    amrex::Abort("ABL with anelastic is not fully functional yet");
}

void ABLAnelastic::post_init_actions()
{
    if (!m_is_anelastic) {
        return;
    }
    initialize_data();
}

void ABLAnelastic::post_regrid_actions()
{
    if (!m_is_anelastic) {
        return;
    }
    initialize_data();
}

void ABLAnelastic::initialize_data()
{
    m_density.resize(m_axis, m_mesh.Geom());
    m_pressure.resize(m_axis, m_mesh.Geom());
}
} // namespace amr_wind
