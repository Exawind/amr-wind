#include <AMReX_BC_TYPES.H>
#include "amr-wind/incflo.H"
#include <cmath>

#include "amr-wind/core/Physics.H"
#include "amr-wind/wind_energy/ABL.H"
#include "amr-wind/physics/BoussinesqBubble.H"
#include "amr-wind/utilities/tagging/RefinementCriteria.H"
#include "amr-wind/utilities/tagging/CartBoxRefinement.H"

using namespace amrex;

/** Parse the input file and populate parameters
 */
void incflo::ReadParameters()
{

    { // Prefix incflo
        ParmParse pp("incflo");

        pp.query("verbose", m_verbose);

        pp.query("initial_iterations", m_initial_iterations);
        pp.query("do_initial_proj", m_do_initial_proj);

        // Physics
        pp.query("constant_density", m_constant_density);

        // Godunov-related flags
        pp.query("use_godunov", m_use_godunov);

        // The default for diffusion_type is 2, i.e. the default m_diff_type is
        // DiffusionType::Implicit
        int diffusion_type = 1;
        pp.query("diffusion_type", diffusion_type);
        if (diffusion_type == 0) {
            m_diff_type = DiffusionType::Explicit;
        } else if (diffusion_type == 1) {
            m_diff_type = DiffusionType::Crank_Nicolson;
        } else if (diffusion_type == 2) {
            m_diff_type = DiffusionType::Implicit;
        } else {
            amrex::Abort(
                "We currently require diffusion_type = 0 for explicit, 1 for "
                "Crank-Nicolson or 2 for implicit");
        }

        if (!m_use_godunov && m_time.max_cfl() > 0.5) {
            amrex::Abort(
                "We currently require cfl <= 0.5 when using the MOL advection "
                "scheme");
        }
        if (m_use_godunov && m_time.max_cfl() > 1.0) {
            amrex::Abort(
                "We currently require cfl <= 1.0 when using the Godunov "
                "advection scheme");
        }

    } // end prefix incflo
}

/** Perform initial pressure iterations
 *
 *  Performs a user-defined number of iterations to compute the pressure based
 *  on the initial conditions. This method is only invoked for new simulations
 *  and skipped for restarted simulations from a checkpoint file.
 */
void incflo::InitialIterations()
{
    BL_PROFILE("amr-wind::incflo::InitialIterations()");
    amrex::Print() << "Begin initial pressure iterations. Num. iters = "
                   << m_initial_iterations << std::endl;

    bool explicit_diffusion = (m_diff_type == DiffusionType::Explicit);
    ComputeDt(explicit_diffusion);

    {
        auto& vel = icns().fields().field;
        vel.copy_state(amr_wind::FieldState::Old, amr_wind::FieldState::New);
        vel.state(amr_wind::FieldState::Old).fillpatch(m_time.current_time());

        if (m_constant_density) {
            auto& rho = density();
            rho.copy_state(
                amr_wind::FieldState::Old, amr_wind::FieldState::New);
            rho.state(amr_wind::FieldState::Old)
                .fillpatch(m_time.current_time());
        }

        for (auto& eqn : scalar_eqns()) {
            auto& scal = eqn->fields().field;
            scal.copy_state(
                amr_wind::FieldState::Old, amr_wind::FieldState::New);
            scal.state(amr_wind::FieldState::Old)
                .fillpatch(m_time.current_time());
        }
    }

    for (int iter = 0; iter < m_initial_iterations; ++iter) {
        if (m_verbose != 0) {
            amrex::Print() << "In initial_iterations: iter = " << iter << "\n";
        }

        ApplyPredictor(true);

        {
            auto& vel = icns().fields().field;
            // ensure velocity is in stretched mesh space
            if (vel.in_uniform_space() && m_sim.has_mesh_mapping()) {
                vel.to_stretched_space();
            }
            vel.copy_state(
                amr_wind::FieldState::New, amr_wind::FieldState::Old);

            if (m_constant_density) {
                auto& rho = density();
                rho.copy_state(
                    amr_wind::FieldState::New, amr_wind::FieldState::Old);
            }

            for (auto& eqn : scalar_eqns()) {
                auto& scal = eqn->fields().field;
                scal.copy_state(
                    amr_wind::FieldState::New, amr_wind::FieldState::Old);
            }
        }
    }

    // Add mean pressure back if available
    if (m_repo.field_exists("reference_pressure")) {
        auto& pressure = m_repo.get_field("p");
        auto& p0 = m_repo.get_field("reference_pressure");
        for (int lev = 0; lev <= finest_level; lev++) {
            amrex::MultiFab::Add(
                pressure(lev), p0(lev), 0, 0, 1, pressure.num_grow()[0]);
        }
    }
    amrex::Print() << "Completed initial pressure iterations" << std::endl
                   << std::endl;
}

/** Ensure initial velocity field is divergence-free.
 *
 *  Performs a \ref incflo::ApplyProjection "projection" step to ensure that the
 *  user-provided initial velocity field is divergence free. This method is only
 *  invoked for new simulations. For restarted simulations using a checkpoint
 *  file, this is not necessary.
 */
void incflo::InitialProjection()
{
    BL_PROFILE("amr-wind::incflo::InitialProjection()");

    amrex::Print() << "Begin initial projection" << std::endl;
    if (m_verbose != 0) {
        PrintMaxValues("before initial projection");
    }

    Real dummy_dt = 1.0;
    bool incremental = false;
    ApplyProjection(
        density().vec_const_ptrs(), m_time.current_time(), dummy_dt,
        incremental);

    // We set p and gp back to zero (p0 may still be still non-zero)
    pressure().setVal(0.0);
    grad_p().setVal(0.0);

    if (m_verbose != 0) {
        PrintMaxValues("after initial projection");
    }
    amrex::Print() << "Completed initial projection" << std::endl << std::endl;
}
