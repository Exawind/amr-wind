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
    pp.query("reinit_iterations", m_niterations);
    pp.query("reinit_convg_interval", m_calcconvint);
    pp.query("reinit_convg_tolerance", m_tol);
    pp.query("reinit_rlscale", m_niterations);
    pp.query("reinit_upw_margin", m_margin);
    pp.query("reinit_target_cutoff", m_target_cutoff);

    // Queries for coupling options
    pp.query("replace_gradp_postsolve", m_replace_gp);
    // OversetOps does not control these coupling options, merely reports them
    pp.query("disable_coupled_nodal_proj", m_disable_nodal_proj);
    pp.query("disable_coupled_mac_proj", m_disable_mac_proj);

    // Verbosity
    pp.query("verbose", m_verbose);

    // Check for perturbational pressure
    // (will be removed soon)
    amrex::ParmParse pp_icns("ICNS");
    pp_icns.query("use_perturb_pressure", m_perturb_p);

    // Check for vof to determine if multiphase sim
    m_vof_exists = (*m_sim_ptr).repo().field_exists("vof");

    // Set up pointer to MultiPhase physics
    if (m_vof_exists) {
        m_mphase = &(*m_sim_ptr).physics_manager().get<MultiPhase>();
    }

    // Set up field to store pressure gradient
    if (m_replace_gp) {
        m_gp_copy = &(*m_sim_ptr).repo().declare_field("gp_copy", 3);
    }

    // Output parameters if verbose
    parameter_output();
}

void OversetOps::pre_advance_work()
{
    // Pressure gradient not updated for current multiphase approach
    if (!(m_vof_exists && m_use_hs_pgrad)) {
        // Update pressure gradient using updated overset pressure field
        update_gradp();
    }

    if (m_vof_exists) {
        // Reinitialize fields
        sharpen_nalu_data();
        if (m_use_hs_pgrad) {
            // Use hydrostatic pressure gradient
            set_hydrostatic_gradp();
        }
    }

    // If pressure gradient will be replaced, store current pressure gradient
    if (m_replace_gp) {
        auto& gp = (*m_sim_ptr).repo().get_field("gp");
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
#pragma omp parallel if (Gpu::notInLaunchRegion())
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
    if (m_replace_gp) {
        replace_masked_gradp();
    }
}

/* ----------------------------------------------- */
/* PUBLIC FUNCTIONS ABOVE, PRIVATE FUNCTIONS BELOW */
/* ----------------------------------------------- */

void OversetOps::parameter_output()
{
    // Print the details
    if (m_verbose > 0) {
        // Important parameters
        amrex::Print() << "Overset Coupling Parameters: \n"
                       << "---- Coupled nodal projection : "
                       << !m_disable_nodal_proj << std::endl
                       << "---- Coupled MAC projection   : "
                       << !m_disable_mac_proj << std::endl
                       << "---- Replace overset pres grad: " << m_replace_gp
                       << std::endl;
        if (m_vof_exists) {
            amrex::Print() << "Overset Reinitialization Parameters:\n"
                           << "---- Maximum iterations   : " << m_niterations
                           << std::endl
                           << "---- Convergence tolerance: " << m_tol
                           << std::endl
                           << "---- Relative length scale: " << m_rlscale
                           << std::endl
                           << "---- Upwinding VOF margin : " << m_margin
                           << std::endl;
        }
    }
    if (m_verbose > 1 && m_vof_exists) {
        // Less important or less used parameters
        amrex::Print() << "---- Calc. conv. interval : " << m_calcconvint
                       << std::endl
                       << "---- Target field cutoff  : " << m_target_cutoff
                       << std::endl;
    }
}

void OversetOps::sharpen_nalu_data()
{
    auto& repo = (*m_sim_ptr).repo();
    auto nlevels = repo.num_active_levels();
    auto geom = (*m_sim_ptr).mesh().Geom();

    // Get blanking for cells
    auto& iblank_cell = repo.get_int_field("iblank_cell");

    // Get fields that will be modified
    auto& vof = repo.get_field("vof");
    auto& levelset = repo.get_field("levelset");
    auto& rho = repo.get_field("density");
    auto& velocity = repo.get_field("velocity");

    // Create scratch fields for fluxes
    // 5 components are vof, density, and 3 of velocity
    auto flux_x = repo.create_scratch_field(5, 0, amr_wind::FieldLoc::XFACE);
    auto flux_y = repo.create_scratch_field(5, 0, amr_wind::FieldLoc::YFACE);
    auto flux_z = repo.create_scratch_field(5, 0, amr_wind::FieldLoc::ZFACE);
    // Create scratch field for approximate signed distance function and grad
    // (components 0-2 are gradient, 3 is asdf)
    auto normal_vec = repo.create_scratch_field(3, vof.num_grow()[0] - 1);

    auto target_vof = repo.create_scratch_field(1, vof.num_grow()[0]);

    // Create levelset field
    for (int lev = 0; lev < nlevels; ++lev) {
        // Thickness used here is user parameter, whatever works best
        auto dx = (geom[lev]).CellSizeArray();
        const amrex::Real i_th = m_rlscale * std::cbrt(dx[0] * dx[1] * dx[2]);

        // Populate approximate signed distance function
        overset_ops::populate_psi(levelset(lev), vof(lev), i_th, m_asdf_tiny);
    }

    // Convert levelset to vof to get target_vof
    m_mphase->levelset2vof(iblank_cell, *target_vof);

    // Process target vof for tiny margins from single-phase
    for (int lev = 0; lev < nlevels; ++lev) {
        // A tolerance of 0 should do nothing
        overset_ops::process_vof((*target_vof)(lev), m_target_cutoff);
    }

    // Replace vof with original values in amr domain
    for (int lev = 0; lev < nlevels; ++lev) {
        overset_ops::harmonize_vof(
            (*target_vof)(lev), vof(lev), iblank_cell(lev));
    }

    // Put fluxes in vector for averaging down during iterations
    amrex::Vector<amrex::Array<amrex::MultiFab*, AMREX_SPACEDIM>> fluxes(
        repo.num_active_levels());
    for (int lev = 0; lev < nlevels; ++lev) {
        fluxes[lev][0] = &(*flux_x)(lev);
        fluxes[lev][1] = &(*flux_y)(lev);
        fluxes[lev][2] = &(*flux_z)(lev);
    }

    // Pseudo-time loop
    amrex::Real err = 100.0 * m_tol;
    int n = 0;
    while (n < m_niterations && err > m_tol) {
        // Increment step counter
        ++n;

        // Determine if convergence error is calculated this step
        bool cconv = n % m_calcconvint == 0;
        // Zero error if being calculated this step
        err = cconv ? 0.0 : err;

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
                (*target_vof)(lev), (*normal_vec)(lev), velocity(lev), m_margin,
                m_mphase->rho1(), m_mphase->rho2());

            // Process fluxes
            overset_ops::process_fluxes(
                (*flux_x)(lev), (*flux_y)(lev), (*flux_z)(lev),
                iblank_cell(lev));

            // Measure convergence to determine if loop can stop
            if (cconv) {
                // Update error at specified interval of steps
                const amrex::Real err_lev = overset_ops::measure_convergence(
                    (*flux_x)(lev), (*flux_y)(lev), (*flux_z)(lev));
                err = amrex::max(err, err_lev);
            }
        }

        // Average down fluxes across levels for consistency
        for (int lev = nlevels - 1; lev > 0; --lev) {
            amrex::IntVect rr =
                geom[lev].Domain().size() / geom[lev - 1].Domain().size();
            amrex::average_down_faces(
                GetArrOfConstPtrs(fluxes[lev]), fluxes[lev - 1], rr,
                geom[lev - 1]);
        }

        // Get pseudo dt (dtau)
        for (int lev = 0; lev < nlevels; ++lev) {
            // Compare vof fluxes to vof in source cells
            // Convergence tolerance determines what size of fluxes matter
            const amrex::Real ptfac_lev = overset_ops::calculate_pseudo_dt_flux(
                (*flux_x)(lev), (*flux_y)(lev), (*flux_z)(lev), vof(lev),
                m_tol);
            ptfac = amrex::min(ptfac, ptfac_lev);
        }
        amrex::ParallelDescriptor::ReduceRealMin(ptfac);

        // Conform pseudo dt (dtau) to pseudo CFL
        ptfac = pCFL * ptfac;

        // Apply fluxes
        for (int lev = 0; lev < nlevels; ++lev) {
            overset_ops::apply_fluxes(
                (*flux_x)(lev), (*flux_y)(lev), (*flux_z)(lev), vof(lev),
                rho(lev), velocity(lev), ptfac, m_vof_tol);
        }

        // Fillpatch for ghost cells
        vof.fillpatch((*m_sim_ptr).time().current_time());
        velocity.fillpatch((*m_sim_ptr).time().current_time());

        // Update density (fillpatch built in)
        m_mphase->set_density_via_vof();

        // Ensure that err is same across processors
        if (cconv) {
            amrex::ParallelDescriptor::ReduceRealMax(err);
        }

        amrex::Print() << "sharpen step " << n << " " << err << " " << m_tol
                       << std::endl;
    }

    // Purely for debugging via visualization, should be removed later
    // Currently set up to overwrite the levelset field (not used as time
    // evolves) with the post-sharpening velocity magnitude
    for (int lev = 0; lev < nlevels; ++lev) {
        overset_ops::equate_field(levelset(lev), velocity(lev));
    }
}

void OversetOps::set_hydrostatic_gradp()
{
    auto& repo = (*m_sim_ptr).repo();
    auto nlevels = repo.num_active_levels();
    auto geom = (*m_sim_ptr).mesh().Geom();

    // Get blanking for cells
    auto& iblank_cell = repo.get_int_field("iblank_cell");

    // Get fields that will be modified or used
    Field* rho0{nullptr};
    auto& rho = repo.get_field("density");
    auto& gp = repo.get_field("gp");
    if (m_perturb_p) {
        rho0 = &((*m_sim_ptr).repo().get_field("reference_density"));
    } else {
        // Point to existing field, won't be used
        rho0 = &rho;
    }

    // Replace initial gp with best guess (hydrostatic)
    for (int lev = 0; lev < nlevels; ++lev) {
        overset_ops::replace_gradp_hs(
            gp(lev), rho(lev), (*rho0)(lev), iblank_cell(lev),
            m_mphase->gravity()[2], m_perturb_p);
    }
}

void OversetOps::replace_masked_gradp()
{
    auto& repo = (*m_sim_ptr).repo();
    auto nlevels = repo.num_active_levels();

    // Get timestep
    const amrex::Real dt = (*m_sim_ptr).time().deltaT();

    // Get blanking for cells
    auto& iblank_cell = repo.get_int_field("iblank_cell");

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