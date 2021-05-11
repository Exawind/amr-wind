#include "amr-wind/physics/multiphase/MultiPhase.H"
#include "amr-wind/equation_systems/vof/volume_fractions.H"
#include "amr-wind/CFDSim.H"
#include "AMReX_ParmParse.H"
#include "amr-wind/fvm/filter.H"
#include "amr-wind/core/field_ops.H"
#include "amr-wind/equation_systems/BCOps.H"
#include <AMReX_MultiFabUtil.H>
#include "amr-wind/core/SimTime.H"

namespace amr_wind {

MultiPhase::MultiPhase(CFDSim& sim)
    : m_sim(sim)
    , m_velocity(sim.pde_manager().icns().fields().field)
    , m_density(sim.repo().get_field("density"))
{
    amrex::ParmParse pp_multiphase("MultiPhase");
    pp_multiphase.query("interface_capturing_method", m_interface_model);
    pp_multiphase.query("density_fluid1", m_rho1);
    pp_multiphase.query("density_fluid2", m_rho2);
    pp_multiphase.query("verbose", m_verbose);

    // Register either the VOF or levelset equation
    if (amrex::toLower(m_interface_model) == "vof") {
        m_interface_capturing_method = amr_wind::InterfaceCapturingMethod::VOF;
        auto& vof_eqn = sim.pde_manager().register_transport_pde("VOF");
        m_vof = &(vof_eqn.fields().field);
        // Create levelset as a auxilliary field only !
        m_levelset = &(m_sim.repo().get_field("levelset"));
        const amrex::Real levelset_default = 0.0;
        BCScalar bc_ls(*m_levelset);
        bc_ls(levelset_default);
        m_levelset->fillpatch(sim.time().current_time());
    } else if (amrex::toLower(m_interface_model) == "levelset") {
        m_interface_capturing_method = amr_wind::InterfaceCapturingMethod::LS;
        auto& levelset_eqn =
            sim.pde_manager().register_transport_pde("Levelset");
        m_levelset = &(levelset_eqn.fields().field);
    } else {
        amrex::Print() << "Please select an interface capturing model between "
                          "VOF and Levelset: defaultin to VOF "
                       << std::endl;
        m_interface_capturing_method = amr_wind::InterfaceCapturingMethod::VOF;
        auto& vof_eqn = sim.pde_manager().register_transport_pde("VOF");
        m_vof = &(vof_eqn.fields().field);
        // Create levelset as a auxilliary field only !
        m_levelset = &(m_sim.repo().get_field("levelset"));
        const amrex::Real levelset_default = 0.0;
        BCScalar bc_ls(*m_levelset);
        bc_ls(levelset_default);
    }
}

InterfaceCapturingMethod MultiPhase::interface_capturing_method()
{
    return m_interface_capturing_method;
}

void MultiPhase::post_init_actions()
{
    switch (m_interface_capturing_method) {
    case InterfaceCapturingMethod::VOF:
        levelset2vof();
        set_density_via_vof();
        break;
    case InterfaceCapturingMethod::LS:
        set_density_via_levelset();
        break;
    };
}

void MultiPhase::post_advance_work()
{
    switch (m_interface_capturing_method) {
    case InterfaceCapturingMethod::VOF:
        set_density_via_vof();
        // Compute the print the total volume fraction
        if (m_verbose > 0) {
            m_total_volfrac = volume_fraction_sum();
            const auto& geom = m_sim.mesh().Geom();
            const amrex::Real total_vol = geom[0].ProbDomain().volume();
            amrex::Print() << "Volume of Fluid diagnostics:" << std::endl;
            amrex::Print() << "   Water Volume Fractions Sum : "
                           << m_total_volfrac << std::endl;
            amrex::Print() << "   Air Volume Fractions Sum : "
                           << total_vol - m_total_volfrac << std::endl;
            amrex::Print() << " " << std::endl;
        }
        break;
    case InterfaceCapturingMethod::LS:
        set_density_via_levelset();
        break;
    };
}

amrex::Real MultiPhase::volume_fraction_sum()
{
    using namespace amrex;
    BL_PROFILE("amr-wind::multiphase::ComputeVolumeFractionSum");
    const int nlevels = m_sim.repo().num_active_levels();
    const auto& geom = m_sim.mesh().Geom();
    auto& mesh = m_sim.mesh();

    amrex::Real total_volume_frac = 0.0;

    for (int lev = 0; lev < nlevels; ++lev) {

        amrex::iMultiFab level_mask;
        if (lev < nlevels - 1) {
            level_mask = makeFineMask(
                mesh.boxArray(lev), mesh.DistributionMap(lev),
                mesh.boxArray(lev + 1), amrex::IntVect(2), 1, 0);
        } else {
            level_mask.define(
                mesh.boxArray(lev), mesh.DistributionMap(lev), 1, 0,
                amrex::MFInfo());
            level_mask.setVal(1);
        }

        auto& vof = (*m_vof)(lev);
        const amrex::Real cell_vol = geom[lev].CellSize()[0] *
                                     geom[lev].CellSize()[1] *
                                     geom[lev].CellSize()[2];

        total_volume_frac += amrex::ReduceSum(
            vof, level_mask, 0,
            [=] AMREX_GPU_HOST_DEVICE(
                amrex::Box const& bx,
                amrex::Array4<amrex::Real const> const& volfrac,
                amrex::Array4<int const> const& mask_arr) -> amrex::Real {
                amrex::Real vol_fab = 0.0;
                amrex::Loop(bx, [=, &vol_fab](int i, int j, int k) noexcept {
                    vol_fab += volfrac(i, j, k) * mask_arr(i, j, k) * cell_vol;
                });
                return vol_fab;
            });
    }
    amrex::ParallelDescriptor::ReduceRealSum(total_volume_frac);

    return total_volume_frac;
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
    m_density.fillpatch(m_sim.time().current_time());
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
    m_density.fillpatch(m_sim.time().current_time());
}

// Reconstructing the volume fraction with the levelset
void MultiPhase::levelset2vof()
{
    const int nlevels = m_sim.repo().num_active_levels();
    (*m_levelset).fillpatch(m_sim.time().new_time());
    const auto& geom = m_sim.mesh().Geom();

    for (int lev = 0; lev < nlevels; ++lev) {
        auto& levelset = (*m_levelset)(lev);
        auto& vof = (*m_vof)(lev);
        const auto& dx = geom[lev].CellSizeArray();

        for (amrex::MFIter mfi(levelset); mfi.isValid(); ++mfi) {
            const auto& vbx = mfi.validbox();
            const amrex::Array4<amrex::Real>& phi = levelset.array(mfi);
            const amrex::Array4<amrex::Real>& volfrac = vof.array(mfi);
            const amrex::Real eps = 2. * std::cbrt(dx[0] * dx[1] * dx[2]);
            amrex::ParallelFor(
                vbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    amrex::Real mx, my, mz;
                    multiphase::youngs_fd_normal(i, j, k, phi, mx, my, mz);
                    mx = std::abs(mx / 32.);
                    my = std::abs(my / 32.);
                    mz = std::abs(mz / 32.);
                    amrex::Real normL1 = (mx + my + mz);
                    mx = mx / normL1;
                    my = my / normL1;
                    mz = mz / normL1;
                    // Make sure that alpha is negative far away from the
                    // interface
                    amrex::Real alpha;
                    if (phi(i, j, k) < -eps) {
                        alpha = -1.0;
                    } else {
                        alpha = phi(i, j, k) / normL1;
                        alpha = alpha + 0.5;
                    }
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
