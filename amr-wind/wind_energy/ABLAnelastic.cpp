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
    auto& rho0 = m_sim.repo().declare_field("reference_density", 1, 0, 1);
    auto& p0 = m_sim.repo().declare_nd_field("reference_pressure", 1, 0, 1);

    rho0.setVal(m_rho0_const);
    p0.setVal(m_atmospheric_pressure);

    const auto axis = m_axis;
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> gravity{
        {m_gravity[0], m_gravity[1], m_gravity[2]}};
    for (int lev = 0; lev < p0.repo().num_active_levels(); ++lev) {
        const auto& darrs = rho0(lev).const_arrays();
        const auto& parrs = p0(lev).arrays();
        const amrex::IntVect ngs(0);
        amrex::ParallelFor(
            rho0(lev), ngs,
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
                const auto& darr = darrs[nbx];
                const auto& parr = parrs[nbx];
                const amrex::IntVect iv(i, j, k);
                const auto ivp = iv + amrex::IntVect::TheDimensionVector(axis);
                parr(ivp) = parr(iv) - darr(iv) * gravity[m_axis];
            });
        amrex::Gpu::synchronize();
    }
}
} // namespace amr_wind
