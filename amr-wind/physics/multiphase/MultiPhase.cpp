#include "amr-wind/physics/multiphase/MultiPhase.H"
#include "amr-wind/equation_systems/vof/VolumeFractions_K.H"
#include "amr-wind/CFDSim.H"
#include "AMReX_ParmParse.H"
#include "amr-wind/fvm/gradient.H"
#include "amr-wind/core/field_ops.H"

namespace amr_wind {

MultiPhase::MultiPhase(CFDSim& sim)
    : m_sim(sim)
    , m_velocity(sim.pde_manager().icns().fields().field)
    , m_density(sim.repo().get_field("density"))
{
    amrex::ParmParse pp_multiphase("MultiPhase");
    pp_multiphase.query("interface_tracking_method", m_interface_model);
    pp_multiphase.query("density_fluid1", m_rho1);
    pp_multiphase.query("density_fluid2", m_rho2);

    // Register either the VOF or levelset equation
    if (amrex::toLower(m_interface_model) == "vof") {
        m_interface_tracking_method = amr_wind::InterfaceTrackingMethod::VOF;
        auto& vof_eqn = sim.pde_manager().register_transport_pde("VOF");
        m_vof = &(vof_eqn.fields().field);
        // Create levelset as a auxilliary field only !
        m_levelset = &(m_sim.repo().get_field("levelset"));
    } else if (amrex::toLower(m_interface_model) == "levelset") {
        m_interface_tracking_method = amr_wind::InterfaceTrackingMethod::LS;
        auto& levelset_eqn =
            sim.pde_manager().register_transport_pde("Levelset");
        m_levelset = &(levelset_eqn.fields().field);
    } else {
        amrex::Print() << "Please select an interface tracking model between "
                          "VOF and Levelset: defaultin to VOF "
                       << std::endl;
        m_interface_tracking_method = amr_wind::InterfaceTrackingMethod::VOF;
        auto& vof_eqn = sim.pde_manager().register_transport_pde("VOF");
        m_vof = &(vof_eqn.fields().field);
        // Create levelset as a auxilliary field only !
        m_levelset = &(m_sim.repo().get_field("levelset"));
    }
}

InterfaceTrackingMethod MultiPhase::interface_tracking_method()
{
    return m_interface_tracking_method;
}

void MultiPhase::post_init_actions()
{

    if (m_interface_tracking_method == InterfaceTrackingMethod::VOF) {
        levelset2vof();
        set_density_via_vof();
    } else if (m_interface_tracking_method == InterfaceTrackingMethod::LS) {
        set_density_via_levelset();
    }
    m_density.fillpatch(m_sim.time().current_time());
}

void MultiPhase::post_advance_work()
{
    if (m_interface_tracking_method == InterfaceTrackingMethod::VOF) {
        set_density_via_vof();
    } else if (m_interface_tracking_method == InterfaceTrackingMethod::LS) {
        set_density_via_levelset();
    }
    m_density.fillpatch(m_sim.time().current_time());
}

void MultiPhase::set_density_via_levelset()
{
    const int nlevels = m_sim.repo().num_active_levels();
    const auto& geom = m_sim.mesh().Geom();

    for (int lev = 0; lev < nlevels; ++lev) {
        auto& density = m_density(lev);
        auto& levelset = (*m_levelset)(lev);

        for (amrex::MFIter mfi(density); mfi.isValid(); ++mfi) {
            const auto& vbx = mfi.validbox();
            const auto& dx = geom[lev].CellSizeArray();

            const amrex::Array4<amrex::Real>& phi = levelset.array(mfi);
            const amrex::Array4<amrex::Real>& rho = density.array(mfi);
            const amrex::Real eps = std::cbrt(2. * dx[0] * dx[1] * dx[2]);
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
                        H = 0.5 *
                            (1.0 + phi(i, j, k) / eps +
                             1.0 / M_PI * std::sin(phi(i, j, k) * M_PI / eps));
                    }
                    rho(i, j, k) = rho1 * H + rho2 * (1.0 - H);
                });
        }
    }
}

void MultiPhase::set_density_via_vof()
{
    const int nlevels = m_sim.repo().num_active_levels();

    for (int lev = 0; lev < nlevels; ++lev) {
        auto& density = m_density(lev);
        auto& vof = (*m_vof)(lev);

        for (amrex::MFIter mfi(density); mfi.isValid(); ++mfi) {
            const auto& vbx = mfi.validbox();
            const amrex::Array4<amrex::Real>& F = vof.array(mfi);
            const amrex::Array4<amrex::Real>& rho = density.array(mfi);
            const amrex::Real rho1 = m_rho1;
            const amrex::Real rho2 = m_rho2;
            amrex::ParallelFor(
                vbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    rho(i, j, k) =
                        rho1 * F(i, j, k) + rho2 * (1.0 - F(i, j, k));
                });
        }
    }
}

// Reconstructing the volume fraction with the levelset
void MultiPhase::levelset2vof()
{
    const int nlevels = m_sim.repo().num_active_levels();
    const auto& geom = m_sim.mesh().Geom();

    auto& normal = m_sim.repo().get_field("interface_normal");
    (*m_levelset).fillpatch(m_sim.time().new_time());
    fvm::gradient(normal, (*m_levelset));
    normal.fillpatch(m_sim.time().new_time());
    field_ops::normalize(normal);

    for (int lev = 0; lev < nlevels; ++lev) {
        auto& levelset = (*m_levelset)(lev);
        auto& vof = (*m_vof)(lev);
        for (amrex::MFIter mfi(levelset); mfi.isValid(); ++mfi) {
            const auto& vbx = mfi.validbox();
            const auto& dx = geom[lev].CellSizeArray();

            const amrex::Array4<amrex::Real>& phi = levelset.array(mfi);
            const amrex::Array4<amrex::Real>& F = vof.array(mfi);
            const amrex::Array4<amrex::Real>& n = normal(lev).array(mfi);

            amrex::ParallelFor(
                vbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    // Do a linear recontruction of the interface based on least
                    // squares
                    n(i, j, k, 0) = -n(i, j, k, 0) * dx[0];
                    n(i, j, k, 1) = -n(i, j, k, 1) * dx[1];
                    n(i, j, k, 2) = -n(i, j, k, 2) * dx[2];
                    amrex::Real mx = amrex::Math::abs(n(i, j, k, 0));
                    amrex::Real my = amrex::Math::abs(n(i, j, k, 1));
                    amrex::Real mz = amrex::Math::abs(n(i, j, k, 2));
                    amrex::Real normL1 = mx + my + mz;
                    mx = mx / normL1;
                    my = my / normL1;
                    mz = mz / normL1;
                    amrex::Real alpha = phi(i, j, k) / normL1;
                    alpha = alpha + 0.5;
                    if (alpha >= 1.0) {
                        F(i, j, k) = 1.0;
                    } else if (alpha <= 0.0) {
                        F(i, j, k) = 0.0;
                    } else {
                        F(i, j, k) = FL3D(mx, my, mz, alpha, 0.0, 1.0);
                    }
                });
        }
    }
}

} // namespace amr_wind
