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
        // Create levelset as a auxiliary field only !
        m_levelset = &(m_sim.repo().get_field("levelset"));
        const amrex::Real levelset_default = 0.0;
        BCFillPatchExtrap bc_ls(*m_levelset);
        bc_ls(levelset_default);
        m_levelset->fillpatch(sim.time().current_time());
    } else if (amrex::toLower(m_interface_model) == "levelset") {
        m_interface_capturing_method = amr_wind::InterfaceCapturingMethod::LS;
        auto& levelset_eqn =
            sim.pde_manager().register_transport_pde("Levelset");
        m_levelset = &(levelset_eqn.fields().field);
    } else {
        amrex::Print() << "Please select an interface capturing model between "
                          "VOF and Levelset: defaulting to VOF "
                       << std::endl;
        m_interface_capturing_method = amr_wind::InterfaceCapturingMethod::VOF;
        auto& vof_eqn = sim.pde_manager().register_transport_pde("VOF");
        m_vof = &(vof_eqn.fields().field);
        // Create levelset as a auxiliary field only !
        m_levelset = &(m_sim.repo().get_field("levelset"));
        const amrex::Real levelset_default = 0.0;
        BCScalar bc_ls(*m_levelset);
        bc_ls(levelset_default);
    }

    // Address pressure approach through input values
    amrex::ParmParse pp_icns("ICNS");
    pp_icns.query("use_perturb_pressure", m_use_perturb_pressure);
    pp_icns.query("reconstruct_true_pressure", m_reconstruct_true_pressure);
    // Declare fields
    if (m_use_perturb_pressure) {
        m_sim.repo().declare_field("reference_density", 1, 0, 1);
        if (m_reconstruct_true_pressure) {
            m_sim.repo().declare_nd_field(
                "reference_pressure", 1, (*m_vof).num_grow()[0], 1);
        }
    }

    // Warn if density specified in single-phase sense
    amrex::ParmParse pp_incflo("incflo");
    if (pp_incflo.contains("density")) {
        amrex::Print() << "WARNING: single-phase density has been specified "
                          "but will not be used! (MultiPhase physics)\n";
    }
    // Always populate gravity
    pp_incflo.queryarr("gravity", m_gravity);
}

InterfaceCapturingMethod MultiPhase::interface_capturing_method() const
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

    m_q0 = momentum_sum(0);
    m_q1 = momentum_sum(1);
    m_q2 = momentum_sum(2);
    m_sumvof0 = volume_fraction_sum();

    // Check if water level is specified (from case definition)
    amrex::ParmParse pp_multiphase("MultiPhase");
    bool is_wlev = pp_multiphase.contains("water_level");
    // Abort if no water level specified
    if (m_use_perturb_pressure && !is_wlev) {
        amrex::Abort(
            "Perturbational pressure requested, but physics case does not "
            "specify water level.");
    }
    if (is_wlev) {
        pp_multiphase.get("water_level", m_water_level0);
    }
    // Make rho0 field if both are specified
    if (m_use_perturb_pressure && is_wlev) {
        // Initialize rho0 field for perturbational density, pressure
        auto& rho0 = m_sim.repo().get_field("reference_density");
        hydrostatic::define_rho0(
            rho0, m_rho1, m_rho2, m_water_level0, m_sim.mesh().Geom());

        // Make p0 field if requested
        if (m_reconstruct_true_pressure) {
            // Initialize p0 field for reconstructing p
            auto& p0 = m_sim.repo().get_field("reference_pressure");
            hydrostatic::define_p0(
                p0, m_rho1, m_rho2, m_water_level0, m_gravity[2],
                m_sim.mesh().Geom());
        }
    }
}

void MultiPhase::post_regrid_actions()
{
    // Reinitialize rho0 if needed
    if (m_use_perturb_pressure) {
        auto& rho0 = m_sim.repo().declare_field("reference_density", 1, 0, 1);
        hydrostatic::define_rho0(
            rho0, m_rho1, m_rho2, m_water_level0, m_sim.mesh().Geom());
        // Reinitialize p0 if needed
        if (m_reconstruct_true_pressure) {
            auto ng = (*m_vof).num_grow();
            auto& p0 = m_sim.repo().declare_nd_field(
                "reference_pressure", 1, ng[0], 1);
            hydrostatic::define_p0(
                p0, m_rho1, m_rho2, m_water_level0, m_gravity[2],
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
            amrex::Real mom_x = momentum_sum(0) - m_q0;
            amrex::Real mom_y = momentum_sum(1) - m_q1;
            amrex::Real mom_z = momentum_sum(2) - m_q2;
            const auto& geom = m_sim.mesh().Geom();
            const amrex::Real total_vol = geom[0].ProbDomain().volume();
            amrex::Print() << "Volume of Fluid diagnostics:" << std::endl;
            amrex::Print() << "   Water Volume Fractions Sum, Difference : "
                           << m_total_volfrac << " "
                           << m_total_volfrac - m_sumvof0 << std::endl;
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
                mesh.boxArray(lev + 1), mesh.refRatio(lev), 1, 0);
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
                mesh.boxArray(lev + 1), mesh.refRatio(lev), 1, 0);
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

        const auto& dx = geom[lev].CellSizeArray();

        const auto& phi_arrs = levelset.const_arrays();
        const auto& rho_arrs = density.arrays();
        const amrex::Real eps = std::cbrt(2. * dx[0] * dx[1] * dx[2]);
        const amrex::Real captured_rho1 = m_rho1;
        const amrex::Real captured_rho2 = m_rho2;
        amrex::ParallelFor(
            density, m_density.num_grow(),
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
                amrex::Real smooth_heaviside;
                if (phi_arrs[nbx](i, j, k) > eps) {
                    smooth_heaviside = 1.0;
                } else if (phi_arrs[nbx](i, j, k) < -eps) {
                    smooth_heaviside = 0.;
                } else {
                    smooth_heaviside =
                        0.5 *
                        (1.0 + phi_arrs[nbx](i, j, k) / eps +
                         1.0 / M_PI *
                             std::sin(phi_arrs[nbx](i, j, k) * M_PI / eps));
                }
                rho_arrs[nbx](i, j, k) =
                    captured_rho1 * smooth_heaviside +
                    captured_rho2 * (1.0 - smooth_heaviside);
            });
    }
    amrex::Gpu::streamSynchronize();
}

void MultiPhase::set_density_via_vof(amr_wind::FieldState fstate)
{
    const int nlevels = m_sim.repo().num_active_levels();

    for (int lev = 0; lev < nlevels; ++lev) {
        auto& density = m_density.state(fstate)(lev);
        auto& vof = (*m_vof).state(fstate)(lev);

        const auto& F_arrs = vof.const_arrays();
        const auto& rho_arrs = density.arrays();
        const amrex::Real captured_rho1 = m_rho1;
        const amrex::Real captured_rho2 = m_rho2;
        amrex::ParallelFor(
            density, m_density.num_grow(),
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
                rho_arrs[nbx](i, j, k) =
                    captured_rho1 * F_arrs[nbx](i, j, k) +
                    captured_rho2 * (1.0 - F_arrs[nbx](i, j, k));
            });
    }
    amrex::Gpu::streamSynchronize();
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
    auto& advalpha_x = m_sim.repo().get_field("advalpha_x");
    auto& advalpha_y = m_sim.repo().get_field("advalpha_y");
    auto& advalpha_z = m_sim.repo().get_field("advalpha_z");

    for (int lev = 0; lev < nlevels; ++lev) {

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        for (amrex::MFIter mfi((*m_vof)(lev)); mfi.isValid(); ++mfi) {
            const auto& bx = mfi.tilebox();
            const auto& bxg1 = amrex::grow(bx, 1);
            const auto& xbx = amrex::surroundingNodes(bx, 0);
            const auto& ybx = amrex::surroundingNodes(bx, 1);
            const auto& zbx = amrex::surroundingNodes(bx, 2);

            const auto& aa_x = advalpha_x(lev).array(mfi);
            const auto& aa_y = advalpha_y(lev).array(mfi);
            const auto& aa_z = advalpha_z(lev).array(mfi);

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
        const auto& vel_arrs = velocity.const_arrays();
        const auto& rho_arrs = density.const_arrays();
        const auto& rhou_arrs = mom.arrays();
        amrex::ParallelFor(
            velocity, amrex::IntVect(1), AMREX_SPACEDIM,
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k, int n) noexcept {
                rhou_arrs[nbx](i, j, k, n) =
                    vel_arrs[nbx](i, j, k, n) * rho_arrs[nbx](i, j, k);
            });
    }
    amrex::Gpu::streamSynchronize();

    // Do the filtering
    fvm::filter((*density_filter), m_density);
    fvm::filter((*momentum_filter), (*momentum));

    for (int lev = 0; lev < nlevels; ++lev) {
        auto& velocity = m_velocity(lev);
        auto& vof = (*m_vof)(lev);
        auto& mom_fil = (*momentum_filter)(lev);
        auto& rho_fil = (*density_filter)(lev);
        const auto& vel_arrs = velocity.arrays();
        const auto& volfrac_arrs = vof.const_arrays();
        const auto& rho_u_f_arrs = mom_fil.const_arrays();
        const auto& rho_f_arrs = rho_fil.const_arrays();
        amrex::ParallelFor(
            velocity, amrex::IntVect(0), AMREX_SPACEDIM,
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k, int n) noexcept {
                if (volfrac_arrs[nbx](i, j, k) <= 0.5) {
                    vel_arrs[nbx](i, j, k, n) = rho_u_f_arrs[nbx](i, j, k, n) /
                                                rho_f_arrs[nbx](i, j, k);
                }
            });
    }
    amrex::Gpu::streamSynchronize();
    m_velocity.fillpatch(m_sim.time().current_time());
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

        const auto& phi_arrs = levelset.const_arrays();
        const auto& volfrac_arrs = vof.arrays();
        const amrex::Real eps = 2. * std::cbrt(dx[0] * dx[1] * dx[2]);
        amrex::ParallelFor(
            levelset,
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
                amrex::Real mx, my, mz;
                multiphase::youngs_finite_difference_normal(
                    i, j, k, phi_arrs[nbx], mx, my, mz);
                mx = std::abs(mx / 32.);
                my = std::abs(my / 32.);
                mz = std::abs(mz / 32.);
                const amrex::Real normL1 = (mx + my + mz);
                mx = mx / normL1;
                my = my / normL1;
                mz = mz / normL1;
                // Make sure that alpha is negative far away from the
                // interface
                const amrex::Real alpha =
                    (phi_arrs[nbx](i, j, k) < -eps)
                        ? -1.0
                        : phi_arrs[nbx](i, j, k) / normL1 + 0.5;
                if (alpha >= 1.0) {
                    volfrac_arrs[nbx](i, j, k) = 1.0;
                } else if (alpha <= 0.0) {
                    volfrac_arrs[nbx](i, j, k) = 0.0;
                } else {
                    volfrac_arrs[nbx](i, j, k) =
                        multiphase::cut_volume(mx, my, mz, alpha, 0.0, 1.0);
                }
            });
    }
    amrex::Gpu::streamSynchronize();

    (*m_vof).fillpatch(m_sim.time().current_time());
}

// Do levelset2vof with iblank neumann and into supplied scratch field
void MultiPhase::levelset2vof(
    const IntField& iblank_cell, ScratchField& vof_scr)
{
    const int nlevels = m_sim.repo().num_active_levels();
    (*m_levelset).fillpatch(m_sim.time().current_time());
    const auto& geom = m_sim.mesh().Geom();

    for (int lev = 0; lev < nlevels; ++lev) {
        auto& levelset = (*m_levelset)(lev);
        auto& vof = vof_scr(lev);
        const auto& dx = geom[lev].CellSizeArray();

        const auto& phi_arrs = levelset.const_arrays();
        const auto& volfrac_arrs = vof.arrays();
        const auto& iblank_arrs = iblank_cell(lev).const_arrays();
        const amrex::Real eps = 2. * std::cbrt(dx[0] * dx[1] * dx[2]);
        amrex::ParallelFor(
            levelset,
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
                // Neumann of levelset across iblank boundaries
                int ibdy =
                    (iblank_arrs[nbx](i, j, k) != iblank_arrs[nbx](i - 1, j, k))
                        ? -1
                        : 0;
                int jbdy =
                    (iblank_arrs[nbx](i, j, k) != iblank_arrs[nbx](i, j - 1, k))
                        ? -1
                        : 0;
                int kbdy =
                    (iblank_arrs[nbx](i, j, k) != iblank_arrs[nbx](i, j, k - 1))
                        ? -1
                        : 0;
                // no cell should be isolated such that -1 and 1 are
                // needed
                ibdy =
                    (iblank_arrs[nbx](i, j, k) != iblank_arrs[nbx](i + 1, j, k))
                        ? +1
                        : ibdy;
                jbdy =
                    (iblank_arrs[nbx](i, j, k) != iblank_arrs[nbx](i, j + 1, k))
                        ? +1
                        : jbdy;
                kbdy =
                    (iblank_arrs[nbx](i, j, k) != iblank_arrs[nbx](i, j, k + 1))
                        ? +1
                        : kbdy;
                amrex::Real mx, my, mz;
                multiphase::youngs_finite_difference_normal_neumann(
                    i, j, k, ibdy, jbdy, kbdy, phi_arrs[nbx], mx, my, mz);
                mx = std::abs(mx / 32.);
                my = std::abs(my / 32.);
                mz = std::abs(mz / 32.);
                const amrex::Real normL1 = (mx + my + mz);
                mx = mx / normL1;
                my = my / normL1;
                mz = mz / normL1;
                // Make sure that alpha is negative far away from the
                // interface
                const amrex::Real alpha =
                    (phi_arrs[nbx](i, j, k) < -eps)
                        ? -1.0
                        : phi_arrs[nbx](i, j, k) / normL1 + 0.5;
                if (alpha >= 1.0) {
                    volfrac_arrs[nbx](i, j, k) = 1.0;
                } else if (alpha <= 0.0) {
                    volfrac_arrs[nbx](i, j, k) = 0.0;
                } else {
                    volfrac_arrs[nbx](i, j, k) =
                        multiphase::cut_volume(mx, my, mz, alpha, 0.0, 1.0);
                }
            });
    }
    amrex::Gpu::streamSynchronize();

    vof_scr.fillpatch(m_sim.time().current_time());
}

} // namespace amr_wind
