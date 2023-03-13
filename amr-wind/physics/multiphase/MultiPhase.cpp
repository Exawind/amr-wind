#include "amr-wind/physics/multiphase/MultiPhase.H"
#include "amr-wind/equation_systems/vof/volume_fractions.H"
#include "amr-wind/physics/multiphase/hydrostatic_ops.H"
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
    pp_multiphase.query("interface_smoothing", m_interface_smoothing);
    pp_multiphase.query("interface_smoothing_frequency", m_smooth_freq);

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

    // Address pressure approach through input values
    amrex::ParmParse pp_icns("ICNS");
    pp_icns.query("use_perturb_pressure", is_pptb);
    pp_icns.query("reconstruct_true_pressure", is_ptrue);
    // Declare fields
    if (is_pptb) {
        m_sim.repo().declare_field("reference_density", 1, 0, 1);
        if (is_ptrue) {
            m_sim.repo().declare_nd_field(
                "reference_pressure", 1, (*m_vof).num_grow()[0], 1);
        }
    }
}

InterfaceCapturingMethod MultiPhase::interface_capturing_method()
{
    return m_interface_capturing_method;
}

void MultiPhase::post_init_actions()
{

    const auto& io_mgr = m_sim.io_manager();
    if (!io_mgr.is_restart()) {
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

    q0 = momentum_sum(0);
    q1 = momentum_sum(1);
    q2 = momentum_sum(2);
    sumvof0 = volume_fraction_sum();

    // Check if water level is specified (from case definition)
    amrex::ParmParse pp_multiphase("MultiPhase");
    bool is_wlev = pp_multiphase.contains("water_level");
    // Abort if no water level specified
    if (is_pptb && !is_wlev) {
        amrex::Abort(
            "Perturbational pressure requested, but physics case does not "
            "specify water level.");
    }
    // Make rho0 field if both are specified
    if (is_pptb && is_wlev) {
        pp_multiphase.get("water_level", water_level0);
        // Initialize rho0 field for perturbational density, pressure
        auto& rho0 = m_sim.repo().get_field("reference_density");
        hydrostatic::define_rho0(
            rho0, m_rho1, m_rho2, water_level0, m_sim.mesh().Geom());

        // Make p0 field if requested
        if (is_ptrue) {
            // Initialize p0 field for reconstructing p
            amrex::ParmParse pp("incflo");
            pp.queryarr("gravity", m_gravity);
            auto& p0 = m_sim.repo().get_field("reference_pressure");
            hydrostatic::define_p0(
                p0, m_rho1, m_rho2, water_level0, m_gravity[2],
                m_sim.mesh().Geom());
        }
    }
}

void MultiPhase::post_regrid_actions()
{
    // Reinitialize rho0 if needed
    if (is_pptb) {
        auto& rho0 = m_sim.repo().declare_field("reference_density", 1, 0, 1);
        hydrostatic::define_rho0(
            rho0, m_rho1, m_rho2, water_level0, m_sim.mesh().Geom());
        // Reinitialize p0 if needed
        if (is_ptrue) {
            auto ng = (*m_vof).num_grow();
            auto& p0 = m_sim.repo().declare_nd_field(
                "reference_pressure", 1, ng[0], 1);
            hydrostatic::define_p0(
                p0, m_rho1, m_rho2, water_level0, m_gravity[2],
                m_sim.mesh().Geom());
        }
    }
}

void MultiPhase::pre_advance_work()
{
    switch (m_interface_capturing_method) {
    case InterfaceCapturingMethod::VOF:
        if (m_interface_smoothing &&
            m_sim.time().time_index() % m_smooth_freq == 0) {
            amrex::Print() << "Smoothing the air-sea interface : "
                           << m_sim.time().current_time() << std::endl;
            favre_filtering();
        }
        break;
    case InterfaceCapturingMethod::LS:
        break;
    };
}

void MultiPhase::post_advance_work()
{
    switch (m_interface_capturing_method) {
    case InterfaceCapturingMethod::VOF:
        // Compute and print the total volume fraction, momenta, and differences
        if (m_verbose > 0) {
            m_total_volfrac = volume_fraction_sum();
            amrex::Real mom_x = momentum_sum(0) - q0;
            amrex::Real mom_y = momentum_sum(1) - q1;
            amrex::Real mom_z = momentum_sum(2) - q2;
            const auto& geom = m_sim.mesh().Geom();
            const amrex::Real total_vol = geom[0].ProbDomain().volume();
            amrex::Print() << "Volume of Fluid diagnostics:" << std::endl;
            amrex::Print() << "   Water Volume Fractions Sum, Difference : "
                           << m_total_volfrac << " "
                           << m_total_volfrac - sumvof0 << std::endl;
            amrex::Print() << "   Air Volume Fractions Sum : "
                           << total_vol - m_total_volfrac << std::endl;
            amrex::Print() << "   Total Momentum Difference (x, y, z) : "
                           << mom_x << " " << mom_y << " " << mom_z
                           << std::endl;
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
    const auto& mesh = m_sim.mesh();

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

amrex::Real MultiPhase::momentum_sum(int n)
{
    using namespace amrex;
    BL_PROFILE("amr-wind::multiphase::ComputeVolumeFractionSum");
    const int nlevels = m_sim.repo().num_active_levels();
    const auto& geom = m_sim.mesh().Geom();
    const auto& mesh = m_sim.mesh();

    amrex::Real total_momentum = 0.0;

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

        auto& velocity = m_sim.repo().get_field("velocity")(lev);
        auto& density = m_sim.repo().get_field("density")(lev);
        const amrex::Real cell_vol = geom[lev].CellSize()[0] *
                                     geom[lev].CellSize()[1] *
                                     geom[lev].CellSize()[2];

        total_momentum += amrex::ReduceSum(
            velocity, density, level_mask, 0,
            [=] AMREX_GPU_HOST_DEVICE(
                amrex::Box const& bx,
                amrex::Array4<amrex::Real const> const& vel,
                amrex::Array4<amrex::Real const> const& dens,
                amrex::Array4<int const> const& mask_arr) -> amrex::Real {
                amrex::Real vol_fab = 0.0;
                amrex::Loop(bx, [=, &vol_fab](int i, int j, int k) noexcept {
                    vol_fab += vel(i, j, k, n) * dens(i, j, k) *
                               mask_arr(i, j, k) * cell_vol;
                });
                return vol_fab;
            });
    }
    amrex::ParallelDescriptor::ReduceRealSum(total_momentum);

    return total_momentum;
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
            const amrex::Real captured_rho1 = m_rho1;
            const amrex::Real captured_rho2 = m_rho2;
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
                    rho(i, j, k) = captured_rho1 * smooth_heaviside +
                                   captured_rho2 * (1.0 - smooth_heaviside);
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
            const amrex::Real captured_rho1 = m_rho1;
            const amrex::Real captured_rho2 = m_rho2;
            amrex::ParallelFor(
                vbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    rho(i, j, k) = captured_rho1 * F(i, j, k) +
                                   captured_rho2 * (1.0 - F(i, j, k));
                });
        }
    }
    m_density.fillpatch(m_sim.time().current_time());
}

void MultiPhase::set_nph_density()
{

    amr_wind::field_ops::lincomb(
        m_density.state(amr_wind::FieldState::NPH), 0.5,
        m_density.state(amr_wind::FieldState::Old), 0, 0.5, m_density, 0, 0,
        m_density.num_comp(), 1);
}

// Using phase densities to convert advected VOF arrays to advected density
void MultiPhase::calculate_advected_facedensity()
{
    const int nlevels = m_sim.repo().num_active_levels();
    amrex::Real c_r1 = m_rho1;
    amrex::Real c_r2 = m_rho2;

    // Get advected vof terms at each face
    // cppcheck-suppress constVariable
    auto& advalpha_x = m_sim.repo().get_field("advalpha_x");
    // cppcheck-suppress constVariable
    auto& advalpha_y = m_sim.repo().get_field("advalpha_y");
    // cppcheck-suppress constVariable
    auto& advalpha_z = m_sim.repo().get_field("advalpha_z");

    for (int lev = 0; lev < nlevels; ++lev) {

        for (amrex::MFIter mfi((*m_vof)(lev)); mfi.isValid(); ++mfi) {
            const auto& bx = mfi.tilebox();
            const auto& bxg1 = amrex::grow(bx, 1);
            const auto& xbx = amrex::surroundingNodes(bx, 0);
            const auto& ybx = amrex::surroundingNodes(bx, 1);
            const auto& zbx = amrex::surroundingNodes(bx, 2);

            auto aa_x = advalpha_x(lev).array(mfi);
            auto aa_y = advalpha_y(lev).array(mfi);
            auto aa_z = advalpha_z(lev).array(mfi);

            amrex::ParallelFor(
                bxg1, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    // Volume terms at each face become density terms
                    if (xbx.contains(i, j, k)) {
                        aa_x(i, j, k) =
                            c_r1 * aa_x(i, j, k) + c_r2 * (1.0 - aa_x(i, j, k));
                    }
                    if (ybx.contains(i, j, k)) {
                        aa_y(i, j, k) =
                            c_r1 * aa_y(i, j, k) + c_r2 * (1.0 - aa_y(i, j, k));
                    }
                    if (zbx.contains(i, j, k)) {
                        aa_z(i, j, k) =
                            c_r1 * aa_z(i, j, k) + c_r2 * (1.0 - aa_z(i, j, k));
                    }
                });
        }
    }
}

void MultiPhase::favre_filtering()
{
    const int nlevels = m_sim.repo().num_active_levels();

    // create scratch fields
    auto density_filter =
        m_sim.repo().create_scratch_field(1, 1, FieldLoc::CELL);
    auto momentum =
        m_sim.repo().create_scratch_field(AMREX_SPACEDIM, 1, FieldLoc::CELL);
    auto momentum_filter =
        m_sim.repo().create_scratch_field(AMREX_SPACEDIM, 1, FieldLoc::CELL);

    for (int lev = 0; lev < nlevels; ++lev) {
        auto& mom = (*momentum)(lev);
        auto& velocity = m_velocity(lev);
        auto& density = m_density(lev);

        for (amrex::MFIter mfi(velocity); mfi.isValid(); ++mfi) {
            const auto& bx = mfi.growntilebox(1);
            const amrex::Array4<amrex::Real>& vel = velocity.array(mfi);
            const amrex::Array4<amrex::Real>& rho = density.array(mfi);
            const amrex::Array4<amrex::Real>& rhou = mom.array(mfi);
            amrex::ParallelFor(
                bx, AMREX_SPACEDIM,
                [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
                    rhou(i, j, k, n) = vel(i, j, k, n) * rho(i, j, k);
                });
        }
    }
    // Do the filtering
    fvm::filter((*density_filter), m_density);
    fvm::filter((*momentum_filter), (*momentum));

    for (int lev = 0; lev < nlevels; ++lev) {
        auto& velocity = m_velocity(lev);
        auto& vof = (*m_vof)(lev);
        auto& mom_fil = (*momentum_filter)(lev);
        auto& rho_fil = (*density_filter)(lev);
        for (amrex::MFIter mfi(velocity); mfi.isValid(); ++mfi) {
            const auto& vbx = mfi.validbox();
            const amrex::Array4<amrex::Real>& vel = velocity.array(mfi);
            const amrex::Array4<amrex::Real>& volfrac = vof.array(mfi);
            const amrex::Array4<amrex::Real>& rho_u_f = mom_fil.array(mfi);
            const amrex::Array4<amrex::Real>& rho_f = rho_fil.array(mfi);
            amrex::ParallelFor(
                vbx, AMREX_SPACEDIM,
                [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
                    if (volfrac(i, j, k) <= 0.5) {
                        vel(i, j, k, n) = rho_u_f(i, j, k, n) / rho_f(i, j, k);
                    }
                });
        }
    }
    m_velocity.fillpatch(0.0);
}

// Reconstructing the volume fraction from a levelset function
void MultiPhase::levelset2vof()
{
    const int nlevels = m_sim.repo().num_active_levels();
    (*m_levelset).fillpatch(m_sim.time().current_time());
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
    // Fill ghost and boundary cells before simulation begins
    (*m_vof).fillpatch(0.0);
}

} // namespace amr_wind
