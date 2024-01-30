#include "amr-wind/wind_energy/ABL.H"

namespace amr_wind {

ABLAnelastic::ABLAnelastic(CFDSim& sim) : m_sim(sim)
{
    {
        amrex::ParmParse pp("ABL");
        pp.query("anelastic", m_is_anelastic);
        pp.query("bottom_reference_pressure", m_bottom_reference_pressure);
    }
    std::string godunov_type;
    int conserv = 1;
    {
        amrex::ParmParse pp("incflo");
        pp.queryarr("gravity", m_gravity);
        pp.query("density", m_rho0_const);
        pp.query("godunov_type", godunov_type);
        pp.query("icns_conserv", conserv);
    }
    if (m_is_anelastic) {
        {
            // ensure use of perturbational pressure form
            amrex::ParmParse pp("ICNS");
            pp.add("use_perturb_pressure", (bool)true);
        }
        if (conserv != 1) {
            amrex::Abort("ABLAnelastic is not supported for icns_conserv != 1");
        }
        if (amrex::toLower(godunov_type) == "bds") {
            amrex::Abort("ABLAnelastic is not supported by BDS");
        }
        const auto& density = m_sim.repo().get_field("density");
        auto& ref_density = m_sim.repo().declare_field(
            "reference_density", 1, density.num_grow()[0], 1);
        ref_density.set_default_fillpatch_bc(m_sim.time());
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
    initialize_data();
}

void ABLAnelastic::initialize_data()
{

    m_density.resize(m_axis, m_sim.mesh().Geom());
    m_pressure.resize(m_axis, m_sim.mesh().Geom());

    AMREX_ALWAYS_ASSERT(m_sim.repo().num_active_levels() == m_density.size());
    AMREX_ALWAYS_ASSERT(m_sim.repo().num_active_levels() == m_pressure.size());

    for (int lev = 0; lev < m_density.size(); lev++) {
        auto& dens = m_density.host_data(lev);
        auto& pres = m_pressure.host_data(lev);
        const auto& dx = m_sim.mesh().Geom(lev).CellSizeArray();
        dens.assign(dens.size(), m_rho0_const);
        pres[0] = m_bottom_reference_pressure;
        for (int k = 0; k < pres.size() - 1; k++) {
            pres[k + 1] = pres[k] - dens[k] * m_gravity[m_axis] * dx[m_axis];
        }
    }
    m_density.copy_host_to_device();
    m_pressure.copy_host_to_device();

    auto& rho0 = m_sim.repo().get_field("reference_density");
    auto& p0 = m_sim.repo().get_field("reference_pressure");
    m_density.copy_to_field(rho0);
    m_pressure.copy_to_field(p0);

    rho0.fillpatch(m_sim.time().current_time());
}
} // namespace amr_wind
