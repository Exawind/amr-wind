#include "amr-wind/physics/multiphase/MultiPhase.H"
#include "amr-wind/equation_systems/vof/volume_fractions.H"
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

    // Compute the total volume fraction
    if (m_interface_tracking_method == InterfaceTrackingMethod::VOF) {
        m_total_volfrac = volume_fraction_sum();
        amrex::Print() << "Volume of Fluid 1: " << m_total_volfrac << std::endl;
    }
}

amrex::Real MultiPhase::volume_fraction_sum()
{
    BL_PROFILE("amr-wind::multiphase::ComputeVolumeFractionSum");
    const int nlevels = m_sim.repo().num_active_levels();
    const auto& geom = m_sim.mesh().Geom();

    amrex::Real TotalVolumeFrac = 0.0;

    for (int lev = 0; lev < nlevels; ++lev) {
        auto& vof = (*m_vof)(lev);
        const amrex::Real cell_vol = geom[lev].CellSize()[0] *
                                     geom[lev].CellSize()[1] *
                                     geom[lev].CellSize()[2];

        TotalVolumeFrac += amrex::ReduceSum(
            vof, 0,
            [=] AMREX_GPU_HOST_DEVICE(
                amrex::Box const& bx,
                amrex::Array4<amrex::Real const> const& volfrac)
                -> amrex::Real {
                amrex::Real Vol_Fab = 0.0;
                amrex::Loop(bx, [=, &Vol_Fab](int i, int j, int k) noexcept {
                    Vol_Fab += volfrac(i, j, k) * cell_vol;
                });
                return Vol_Fab;
            });
    }
    amrex::ParallelDescriptor::ReduceRealSum(TotalVolumeFrac);

    return TotalVolumeFrac;
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
                    amrex::Real smooth_heaviside;
                    if (phi(i, j, k) > eps) {
                        smooth_heaviside = 1.0;
                    } else if (phi(i, j, k) < -eps) {
                        smooth_heaviside = 0.;
                    } else {
                        smooth_heaviside =
                            0.5 *
                            (1.0 + phi(i, j, k) / eps +
                             1.0 / M_PI * std::sin(phi(i, j, k) * M_PI / eps));
                    }
                    rho(i, j, k) = rho1 * smooth_heaviside +
                                   rho2 * (1.0 - smooth_heaviside);
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
            const amrex::Array4<amrex::Real>& volfrac = vof.array(mfi);
            const amrex::Array4<amrex::Real>& ifacenorm =
                normal(lev).array(mfi);

            amrex::ParallelFor(
                vbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    // Do a linear recontruction of the interface based on least
                    // squares
                    ifacenorm(i, j, k, 0) = -ifacenorm(i, j, k, 0) * dx[0];
                    ifacenorm(i, j, k, 1) = -ifacenorm(i, j, k, 1) * dx[1];
                    ifacenorm(i, j, k, 2) = -ifacenorm(i, j, k, 2) * dx[2];
                    amrex::Real mx = amrex::Math::abs(ifacenorm(i, j, k, 0));
                    amrex::Real my = amrex::Math::abs(ifacenorm(i, j, k, 1));
                    amrex::Real mz = amrex::Math::abs(ifacenorm(i, j, k, 2));
                    amrex::Real normL1 = mx + my + mz;
                    mx = mx / normL1;
                    my = my / normL1;
                    mz = mz / normL1;
                    amrex::Real alpha = phi(i, j, k) / normL1;
                    alpha = alpha + 0.5;
                    if (alpha >= 1.0) {
                        volfrac(i, j, k) = 1.0;
                    } else if (alpha <= 0.0) {
                        volfrac(i, j, k) = 0.0;
                    } else {
                        volfrac(i, j, k) =
                            multiphase::cut_volume(mx, my, mz, alpha, 0.0, 1.0);
                    }
                });
        }
    }
}

} // namespace amr_wind
