#include "amr-wind/wind_energy/ABL.H"

namespace amr_wind {

ABLAnelastic::ABLAnelastic(CFDSim& sim) : m_mesh(sim.mesh())
{
    {
        amrex::ParmParse pp("ABL");
        pp.query("anelastic", m_is_anelastic);
    }
    {
        amrex::ParmParse pp("incflo");
        pp.queryarr("gravity", m_gravity);
        pp.query("density", m_rho0_const);
    }
    if (m_is_anelastic) {
        amrex::Abort("ABL with anelastic is not fully functional yet");
    }
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

    for (int lev = 0; lev < m_density.size(); lev++) {
        auto dens = m_density.host_data(lev);
        auto pres = m_pressure.host_data(lev);
        dens.assign(dens.size(), m_rho0_const);
        pres[0] = m_atmospheric_pressure;
        for (int k = 0; k < dens.size() - 1; k++) {
            pres[k + 1] =
                pres[k] - 0.5 * (dens[k] + dens[k + 1]) * m_gravity[m_axis];
        }
    }
    m_density.copy_host_to_device();
    m_pressure.copy_host_to_device();
}
} // namespace amr_wind
