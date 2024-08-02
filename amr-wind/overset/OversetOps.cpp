#include "amr-wind/overset/OversetOps.H"
#include "amr-wind/overset/overset_ops_routines.H"
#include "amr-wind/core/field_ops.H"
#include "AMReX_ParmParse.H"
#include "AMReX_MultiFabUtil.H"
#include "amr-wind/core/MLMGOptions.H"
#include "amr-wind/projection/nodal_projection_ops.H"
#include <hydro_NodalProjector.H>

namespace amr_wind {

void OversetOps::initialize(CFDSim& sim)
{
    m_sim_ptr = &sim;
    // Queries for reinitialization options
    amrex::ParmParse pp("Overset");
    pp.query("reinit_iterations", m_n_iterations);
    pp.query("reinit_convg_interval", m_calc_convg_interval);
    pp.query("reinit_convg_tolerance", m_convg_tol);
    pp.query("reinit_rlscale", m_relative_length_scale);
    pp.query("reinit_upw_margin", m_upw_margin);
    pp.query("reinit_target_cutoff", m_target_cutoff);

    // Queries for coupling options
    pp.query("use_hydrostatic_gradp", m_use_hydrostatic_gradp);
    pp.query("replace_gradp_postsolve", m_replace_gradp_postsolve);
    // OversetOps does not control these coupling options, merely reports them
    pp.query("disable_coupled_nodal_proj", m_disable_nodal_proj);
    pp.query("disable_coupled_mac_proj", m_disable_mac_proj);

    pp.query("verbose", m_verbose);

    m_vof_exists = (*m_sim_ptr).repo().field_exists("vof");
    if (m_vof_exists) {
        m_mphase = &(*m_sim_ptr).physics_manager().get<MultiPhase>();

        // Check combination of pressure options
        if (m_mphase->perturb_pressure() &&
            !m_mphase->reconstruct_true_pressure()) {
            amrex::Abort(
                "OversetOps: perturbational pressure is turned on, but true "
                "pressure reconstruction is turned off. This approach will be "
                "incorrect when coupling with Nalu-Wind.");
        }
        if (!m_mphase->perturb_pressure() &&
            m_mphase->reconstruct_true_pressure()) {
            amrex::Print()
                << "WARNING (OversetOps): true pressure reconstruction is "
                   "turned on, but it will remain inactive because "
                   "perturbational pressure is turned off.\n";
        }
    }
    if (m_replace_gradp_postsolve) {
        m_gp_copy = &(*m_sim_ptr).repo().declare_field("gp_copy", 3);
    }

    parameter_output();
}

void OversetOps::pre_advance_work()
{
    if (!(m_vof_exists && m_use_hydrostatic_gradp)) {
        // Update pressure gradient using updated overset pressure field
        update_gradp();
    }

    if (m_vof_exists) {
        // Reinitialize fields
        sharpen_nalu_data();
        if (m_use_hydrostatic_gradp) {
            // Use hydrostatic pressure gradient
            set_hydrostatic_gradp();
        } else {
            if (m_mphase->perturb_pressure()) {
                // Modify to be consistent with internal source terms
                form_perturb_pressure();
            }
            // Update pressure gradient using sharpened pressure field
            update_gradp();
        }
    }

    // If pressure gradient will be replaced, store current pressure gradient
    if (m_replace_gradp_postsolve) {
        const auto& gp = (*m_sim_ptr).repo().get_field("gp");
        for (int lev = 0; lev < (*m_sim_ptr).repo().num_active_levels();
             ++lev) {
            amrex::MultiFab::Copy(
                (*m_gp_copy)(lev), gp(lev), 0, 0, gp(lev).nComp(),
                (m_gp_copy)->num_grow());
        }
    }
}

void OversetOps::update_gradp()
{
    BL_PROFILE("amr-wind::OversetOps::update_gradp");

    // Get relevant fields
    auto& grad_p = (*m_sim_ptr).repo().get_field("gp");
    auto& pressure = (*m_sim_ptr).repo().get_field("p");
    auto& velocity = (*m_sim_ptr).repo().get_field("velocity");

    // Set up projection object
    std::unique_ptr<Hydro::NodalProjector> nodal_projector;
    // Boundary conditions from field, periodicity
    auto bclo = nodal_projection::get_projection_bc(
        amrex::Orientation::low, pressure,
        (*m_sim_ptr).mesh().Geom(0).isPeriodic());
    auto bchi = nodal_projection::get_projection_bc(
        amrex::Orientation::high, pressure,
        (*m_sim_ptr).mesh().Geom(0).isPeriodic());
    // Velocity multifab is needed for proper initialization, but only the size
    // matters for the purpose of calculating gradp, the values do not matter
    int finest_level = (*m_sim_ptr).mesh().finestLevel();
    Vector<MultiFab*> vel;
    for (int lev = 0; lev <= finest_level; ++lev) {
        vel.push_back(&(velocity(lev)));
    }
    amr_wind::MLMGOptions options("nodal_proj");
    // Create nodal projector with unity scaling factor for simplicity
    nodal_projector = std::make_unique<Hydro::NodalProjector>(
        vel, 1.0, (*m_sim_ptr).mesh().Geom(0, finest_level), options.lpinfo());
    // Set MLMG and NodalProjector options
    options(*nodal_projector);
    nodal_projector->setDomainBC(bclo, bchi);

    // Recalculate gradphi with fluxes
    auto gradphi = nodal_projector->calcGradPhi(pressure.vec_ptrs());

    // Transfer pressure gradient to gp field
    for (int lev = 0; lev <= finest_level; lev++) {

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(grad_p(lev), TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            Box const& tbx = mfi.tilebox();
            Array4<Real> const& gp_lev = grad_p(lev).array(mfi);
            Array4<Real const> const& gp_proj = gradphi[lev]->const_array(mfi);
            amrex::ParallelFor(
                tbx, AMREX_SPACEDIM,
                [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
                    gp_lev(i, j, k, n) = gp_proj(i, j, k, n);
                });
        }
    }

    // Averaging down here would be unnecessary; it is built into calcGradPhi
}

void OversetOps::post_advance_work()
{
    // Replace and reapply pressure gradient if requested
    if (m_replace_gradp_postsolve) {
        replace_masked_gradp();
    }
}

/* ----------------------------------------------- */
/* PUBLIC FUNCTIONS ABOVE, PRIVATE FUNCTIONS BELOW */
/* ----------------------------------------------- */

void OversetOps::parameter_output() const
{
    // Print the details
    if (m_verbose > 0) {
        // Important parameters
        amrex::Print() << "\nOverset Coupling Parameters: \n"
                       << "---- Coupled nodal projection : "
                       << !m_disable_nodal_proj << std::endl
                       << "---- Coupled MAC projection   : "
                       << !m_disable_mac_proj << std::endl
                       << "---- Replace overset pres grad: "
                       << m_replace_gradp_postsolve << std::endl;
        if (m_vof_exists) {
            amrex::Print() << "---- Perturbational pressure  : "
                           << m_mphase->perturb_pressure() << std::endl
                           << "---- Reconstruct true pressure: "
                           << m_mphase->reconstruct_true_pressure()
                           << std::endl;
            amrex::Print() << "Overset Reinitialization Parameters:\n"
                           << "---- Maximum iterations   : " << m_n_iterations
                           << std::endl
                           << "---- Convergence tolerance: " << m_convg_tol
                           << std::endl
                           << "---- Relative length scale: "
                           << m_relative_length_scale << std::endl
                           << "---- Upwinding VOF margin : " << m_upw_margin
                           << std::endl;
            if (m_verbose > 1) {
                // Less important or less used parameters
                amrex::Print()
                    << "---- Calc. conv. interval : " << m_calc_convg_interval
                    << std::endl
                    << "---- Target field cutoff  : " << m_target_cutoff
                    << std::endl;
            }
        }
        amrex::Print() << std::endl;
    }
}

void OversetOps::sharpen_nalu_data()
{
    const auto& repo = (*m_sim_ptr).repo();
    const auto nlevels = repo.num_active_levels();
    const auto geom = (*m_sim_ptr).mesh().Geom();

    // Get blanking for cells
    const auto& iblank_cell = repo.get_int_field("iblank_cell");

    // Get fields that will be modified
    auto& vof = repo.get_field("vof");
    auto& levelset = repo.get_field("levelset");
    auto& rho = repo.get_field("density");
    auto& velocity = repo.get_field("velocity");
    auto& gp_noghost = repo.get_field("gp");
    auto& p = repo.get_field("p");

    // 9 components are vof, density, 3 of velocity, 3 of gp, and psource flag
    auto flux_x = repo.create_scratch_field(9, 1, amr_wind::FieldLoc::XFACE);
    auto flux_y = repo.create_scratch_field(9, 1, amr_wind::FieldLoc::YFACE);
    auto flux_z = repo.create_scratch_field(9, 1, amr_wind::FieldLoc::ZFACE);

    auto p_src = repo.create_scratch_field(1, 0, amr_wind::FieldLoc::NODE);
    auto normal_vec = repo.create_scratch_field(3, vof.num_grow()[0] - 1);
    auto target_vof = repo.create_scratch_field(1, vof.num_grow()[0]);

    // Sharpening fluxes (at faces) have 1 ghost, requiring fields to have >= 2
    auto gp_scr = repo.create_scratch_field(3, 2);
    auto& gp = *gp_scr;

    // Give initial max possible value of pseudo-velocity scale
    const auto dx_lev0 = (geom[0]).CellSizeArray();
    const amrex::Real max_pvscale =
        std::min(std::min(dx_lev0[0], dx_lev0[1]), dx_lev0[2]);
    amrex::Real pvscale = max_pvscale;

    // Prep things that do not change with iterations
    for (int lev = 0; lev < nlevels; ++lev) {
        // Thickness used here is user parameter, whatever works best
        const auto dx = (geom[lev]).CellSizeArray();
        const amrex::Real i_th =
            m_relative_length_scale * std::cbrt(dx[0] * dx[1] * dx[2]);

        // Populate approximate signed distance function
        overset_ops::populate_psi(levelset(lev), vof(lev), i_th, m_asdf_tiny);

        gp(lev).setVal(0.0);
        amrex::MultiFab::Copy(gp(lev), gp_noghost(lev), 0, 0, 3, 0);
        gp(lev).FillBoundary(geom[lev].periodicity());

        // Get pseudo-velocity scale, proportional to smallest dx in iblank
        const amrex::Real pvscale_lev =
            overset_ops::calculate_pseudo_velocity_scale(
                iblank_cell(lev), dx, max_pvscale);
        pvscale = std::min(pvscale, pvscale_lev);
    }
    amrex::Gpu::synchronize();
    amrex::ParallelDescriptor::ReduceRealMin(pvscale);

    // Convert levelset to vof to get target_vof
    m_mphase->levelset2vof(iblank_cell, *target_vof);

    // Process target vof for tiny margins from single-phase
    for (int lev = 0; lev < nlevels; ++lev) {
        // A tolerance of 0 should do nothing
        overset_ops::process_vof((*target_vof)(lev), m_target_cutoff);
    }
    amrex::Gpu::synchronize();

    // Replace vof with original values in amr domain
    for (int lev = 0; lev < nlevels; ++lev) {
        overset_ops::harmonize_vof(
            (*target_vof)(lev), vof(lev), iblank_cell(lev));
    }
    amrex::Gpu::synchronize();

    // Put fluxes in vector for averaging down during iterations
    amrex::Vector<amrex::Array<amrex::MultiFab*, AMREX_SPACEDIM>> fluxes(
        repo.num_active_levels());
    for (int lev = 0; lev < nlevels; ++lev) {
        fluxes[lev][0] = &(*flux_x)(lev);
        fluxes[lev][1] = &(*flux_y)(lev);
        fluxes[lev][2] = &(*flux_z)(lev);
    }

    // Pseudo-time loop
    amrex::Real err = 100.0 * m_convg_tol;
    int n = 0;
    while (n < m_n_iterations && err > m_convg_tol) {
        ++n;
        bool calc_convg = n % m_calc_convg_interval == 0;
        // Zero error if being calculated this step
        err = calc_convg ? 0.0 : err;

        // Maximum possible value of pseudo time factor (dtau)
        amrex::Real ptfac = 1.0;
        // Maximum pseudoCFL, 0.5 seems to work well
        const amrex::Real pCFL = 0.5;

        for (int lev = 0; lev < nlevels; ++lev) {
            // Populate normal vector
            overset_ops::populate_normal_vector(
                (*normal_vec)(lev), vof(lev), iblank_cell(lev));

            // Sharpening fluxes for vof, density, and momentum
            overset_ops::populate_sharpen_fluxes(
                (*flux_x)(lev), (*flux_y)(lev), (*flux_z)(lev), vof(lev),
                (*target_vof)(lev), (*normal_vec)(lev), velocity(lev), gp(lev),
                rho(lev), pvscale, m_upw_margin, m_mphase->rho1(),
                m_mphase->rho2());

            // Process fluxes
            overset_ops::process_fluxes_calc_src(
                (*flux_x)(lev), (*flux_y)(lev), (*flux_z)(lev), (*p_src)(lev),
                vof(lev), iblank_cell(lev));

            // Measure convergence to determine if loop can stop
            if (calc_convg) {
                // Update error at specified interval of steps
                const amrex::Real err_lev =
                    overset_ops::measure_convergence(
                        (*flux_x)(lev), (*flux_y)(lev), (*flux_z)(lev)) /
                    pvscale;
                err = amrex::max(err, err_lev);
            }
        }
        amrex::Gpu::synchronize();

        // Average down fluxes across levels for consistency
        for (int lev = nlevels - 1; lev > 0; --lev) {
            amrex::IntVect rr =
                geom[lev].Domain().size() / geom[lev - 1].Domain().size();
            amrex::average_down_faces(
                GetArrOfConstPtrs(fluxes[lev]), fluxes[lev - 1],
                rr, geom[lev - 1]);
        }

        // Get pseudo dt (dtau)
        for (int lev = 0; lev < nlevels; ++lev) {
            const auto dx = (geom[lev]).CellSizeArray();
            // Compare vof fluxes to vof in source cells
            // Convergence tolerance determines what size of fluxes matter
            const amrex::Real ptfac_lev = overset_ops::calculate_pseudo_dt_flux(
                (*flux_x)(lev), (*flux_y)(lev), (*flux_z)(lev), vof(lev), dx,
                m_convg_tol);
            ptfac = amrex::min(ptfac, ptfac_lev);
        }
        amrex::ParallelDescriptor::ReduceRealMin(ptfac);

        // Conform pseudo dt (dtau) to pseudo CFL
        ptfac = pCFL * ptfac;

        // Apply fluxes
        for (int lev = 0; lev < nlevels; ++lev) {
            const auto dx = (geom[lev]).CellSizeArray();

            overset_ops::apply_fluxes(
                (*flux_x)(lev), (*flux_y)(lev), (*flux_z)(lev), (*p_src)(lev),
                vof(lev), rho(lev), velocity(lev), gp(lev), p(lev), dx, ptfac,
                m_vof_tol);

            vof(lev).FillBoundary(geom[lev].periodicity());
            velocity(lev).FillBoundary(geom[lev].periodicity());
            gp(lev).FillBoundary(geom[lev].periodicity());
        }

        amrex::Gpu::synchronize();

        // Update density (fillpatch built in)
        m_mphase->set_density_via_vof();

        // Ensure that err is same across processors
        if (calc_convg) {
            amrex::ParallelDescriptor::ReduceRealMax(err);
        }

        if (m_verbose > 0) {
            amrex::Print() << "OversetOps: sharpen step " << n << "  conv. err "
                           << err << "  tol " << m_convg_tol << std::endl;
        }
    }

    // Fillpatch for pressure to make sure pressure stencil has all points
    p.fillpatch((*m_sim_ptr).time().current_time());

    // Purely for debugging via visualization, should be removed later
    // Currently set up to overwrite the levelset field (not used as time
    // evolves) with the post-sharpening velocity magnitude
    for (int lev = 0; lev < nlevels; ++lev) {
        overset_ops::equate_field(levelset(lev), velocity(lev));
    }
    amrex::Gpu::synchronize();
}

void OversetOps::form_perturb_pressure()
{
    auto& pressure = (*m_sim_ptr).repo().get_field("p");
    const auto& p0 = (*m_sim_ptr).repo().get_field("reference_pressure");
    for (int lev = 0; lev < (*m_sim_ptr).repo().num_active_levels(); lev++) {
        amrex::MultiFab::Subtract(
            pressure(lev), p0(lev), 0, 0, 1, pressure.num_grow()[0]);
    }
}

void OversetOps::set_hydrostatic_gradp()
{
    const auto& repo = (*m_sim_ptr).repo();
    const auto nlevels = repo.num_active_levels();
    const auto geom = (*m_sim_ptr).mesh().Geom();

    // Get blanking for cells
    const auto& iblank_cell = repo.get_int_field("iblank_cell");

    // Get fields that will be modified or used
    Field* rho0{nullptr};
    auto& rho = repo.get_field("density");
    auto& gp = repo.get_field("gp");
    if (m_mphase->perturb_pressure()) {
        rho0 = &((*m_sim_ptr).repo().get_field("reference_density"));
    } else {
        // Point to existing field, won't be used
        rho0 = &rho;
    }

    // Replace initial gp with best guess (hydrostatic)
    for (int lev = 0; lev < nlevels; ++lev) {
        overset_ops::replace_gradp_hydrostatic(
            gp(lev), rho(lev), (*rho0)(lev), iblank_cell(lev),
            m_mphase->gravity()[2], m_mphase->perturb_pressure());
    }
}

void OversetOps::replace_masked_gradp()
{
    const auto& repo = (*m_sim_ptr).repo();
    const auto nlevels = repo.num_active_levels();

    // Get timestep
    const amrex::Real dt = (*m_sim_ptr).time().delta_t();

    // Get blanking for cells
    const auto& iblank_cell = repo.get_int_field("iblank_cell");

    // Get fields that will be modified or used
    auto& vel = repo.get_field("velocity");
    auto& rho = repo.get_field("density");
    auto& gp = repo.get_field("gp");

    // For iblanked cells, replace gp with original gp, to get original vel
    for (int lev = 0; lev < nlevels; ++lev) {
        // Remove pressure gradient term
        overset_ops::apply_pressure_gradient(vel(lev), rho(lev), gp(lev), -dt);
        // Modify pressure gradient
        overset_ops::replace_gradp(
            gp(lev), (*m_gp_copy)(lev), iblank_cell(lev));
        // Reapply pressure gradient term
        overset_ops::apply_pressure_gradient(vel(lev), rho(lev), gp(lev), dt);
    }
}

} // namespace amr_wind
