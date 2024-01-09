#include "amr-wind/wind_energy/ABL.H"

namespace amr_wind {

ABLAnelastic::ABLAnelastic(CFDSim& sim) : m_sim(sim)
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
        m_sim.repo().declare_field("reference_density", 1, 0, 1);
        m_sim.repo().declare_nd_field("reference_pressure", 1, 0, 1);
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
    m_sim.repo().declare_field("reference_density", 1, 0, 1);
    m_sim.repo().declare_nd_field("reference_pressure", 1, 0, 1);
    initialize_data();
}

void ABLAnelastic::initialize_data()
{

    m_density.resize(m_axis, m_sim.mesh().Geom());
    m_pressure.resize(m_axis, m_sim.mesh().Geom());

    for (int lev = 0; lev < m_density.size(); lev++) {
        auto& dens = m_density.host_data(lev);
        auto& pres = m_pressure.host_data(lev);
        dens.assign(dens.size(), m_rho0_const);
        pres[0] = m_atmospheric_pressure;
        const amrex::Real factor =
            lev == 0 ? 1
                     : std::pow(m_sim.mesh().refRatio(lev - 1)[m_axis], lev);
        for (int k = 0; k < pres.size() - 1; k++) {
            pres[k + 1] = pres[k] - dens[k] * m_gravity[m_axis] / factor;
        }
    }
    m_density.copy_host_to_device();
    m_pressure.copy_host_to_device();

    auto& rho0 = m_sim.repo().get_field("reference_density");
    auto& p0 = m_sim.repo().get_field("reference_pressure");
    m_density.to_field(rho0);
    m_pressure.to_field(p0);
}
} // namespace amr_wind
