#include "amr-wind/physics/multiphase/MultiPhase.H"
#include "amr-wind/CFDSim.H"
#include "AMReX_ParmParse.H"

namespace amr_wind {

MultiPhase::MultiPhase(CFDSim& sim)
    : m_sim(sim)
    , m_velocity(sim.pde_manager().icns().fields().field)
    , m_density(sim.repo().get_field("density"))
{
    amrex::ParmParse pp_multiphase("MultiPhase");
    pp_multiphase.query("density_fluid1", m_rho1);
    pp_multiphase.query("density_fluid2", m_rho2);

    // Register both VOF and levelset equations
    auto& vof_eqn = sim.pde_manager().register_transport_pde("VOF");
    auto& levelset_eqn = sim.pde_manager().register_transport_pde("Levelset");
    m_vof = &(vof_eqn.fields().field);
    m_levelset = &(levelset_eqn.fields().field);
}

void MultiPhase::post_init_actions()
{
    const int nlevels = m_sim.repo().num_active_levels();
    const auto& geom = m_sim.mesh().Geom();

    for (int lev = 0; lev < nlevels; ++lev) {
        set_density(lev, geom[lev]);
    }
    m_density.fillpatch(m_sim.time().current_time());
}

void MultiPhase::post_advance_work()
{
    const int nlevels = m_sim.repo().num_active_levels();
    const auto& geom = m_sim.mesh().Geom();
    // Set ve
    for (int lev = 0; lev < nlevels; ++lev) {
        set_density(lev, geom[lev]);
    }
    m_density.fillpatch(m_sim.time().new_time());
}

void MultiPhase::set_density(
    int level, const amrex::Geometry& geom)
{
    auto& density = m_density(level);
    auto& levelset = (*m_levelset)(level);

    for (amrex::MFIter mfi(density); mfi.isValid(); ++mfi) {
        const auto& vbx = mfi.validbox();
        const auto& dx = geom.CellSizeArray();

        const amrex::Array4<amrex::Real>& phi = levelset.array(mfi);
        const amrex::Array4<amrex::Real>& rho = density.array(mfi);
        const amrex::Real eps = std::cbrt(2.*dx[0] * dx[1] * dx[2]);
        const amrex::Real rho1 = m_rho1;
        const amrex::Real rho2 = m_rho2;
        amrex::ParallelFor(
            vbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                amrex::Real H;
                if (phi(i, j, k) > eps) {
                    H = 1.0;
                } else if (phi(i, j, k) < -eps) {
                    H = 0.;
                } else {
                    H = 0.5 * (1 + phi(i, j, k) / (2 * eps) +
                               1. / M_PI * std::sin(phi(i, j, k) * M_PI / eps));
                }
                rho(i, j, k) = rho1 * H + rho2 * (1 - H);
            });
    }
}

} // namespace amr_wind
