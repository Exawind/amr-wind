#include <incflo.H>
#include "Physics.H"
#include <cmath>
#include "field_ops.H"
#include "Godunov.H"
#include "MOL.H"
#include "PDEBase.H"
#include "mac_projection.H"
#include "diffusion.H"
#include "TurbulenceModel.H"
#include "console_io.H"

using namespace amrex;

void incflo::Advance()
{
    BL_PROFILE("amr-wind::incflo::Advance")

    // Compute time step size
    bool explicit_diffusion = (m_diff_type == DiffusionType::Explicit);
    ComputeDt(explicit_diffusion);

    if (m_constant_density) {
        density().advance_states();
        density().state(amr_wind::FieldState::Old).fillpatch(m_time.current_time());
    }

    auto& vel = icns().fields().field;
    vel.advance_states();
    vel.state(amr_wind::FieldState::Old).fillpatch(m_time.current_time());
    for (auto& eqn: scalar_eqns()) {
        auto& field = eqn->fields().field;
        field.advance_states();
        field.state(amr_wind::FieldState::Old).fillpatch(m_time.current_time());
    }

    for (auto& pp: m_sim.physics())
        pp->pre_advance_work();

    ApplyPredictor();

    if (!m_use_godunov) {
        auto& vel = icns().fields().field;
        vel.state(amr_wind::FieldState::New).fillpatch(m_time.current_time());
        for (auto& eqn: scalar_eqns()) {
            auto& field = eqn->fields().field;
            field.state(amr_wind::FieldState::New).fillpatch(m_time.current_time());
        }

        ApplyCorrector();
    }

    if (m_verbose > 1) PrintMaxValues("end of timestep");
}

//
// Apply predictor:
//
//  1. Use u = vel_old to compute
//
//      conv_u  = - u grad u
//      conv_r  = - div( u rho  )
//      conv_t  = - div( u trac )
//      eta_old     = visosity at m_time.current_time()
//      if (m_diff_type == DiffusionType::Explicit)
//         divtau _old = div( eta ( (grad u) + (grad u)^T ) ) / rho^n
//         rhs = u + dt * ( conv + divtau_old )
//      else
//         divtau_old  = 0.0
//         rhs = u + dt * conv
//
//      eta     = eta at new_time
//
//  2. Add explicit forcing term i.e. gravity + lagged pressure gradient
//
//      rhs += dt * ( g - grad(p + p0) / rho^nph )
//
//      Note that in order to add the pressure gradient terms divided by rho,
//      we convert the velocity to momentum before adding and then convert them back.
//
//  3. A. If (m_diff_type == DiffusionType::Implicit)
//        solve implicit diffusion equation for u*
//
//     ( 1 - dt / rho^nph * div ( eta grad ) ) u* = u^n + dt * conv_u
//                                                  + dt * ( g - grad(p + p0) / rho^nph )
//
//     B. If (m_diff_type == DiffusionType::Crank-Nicolson)
//        solve semi-implicit diffusion equation for u*
//
//     ( 1 - (dt/2) / rho^nph * div ( eta_old grad ) ) u* = u^n +
//            dt * conv_u + (dt/2) / rho * div (eta_old grad) u^n
//          + dt * ( g - grad(p + p0) / rho^nph )
//
//  4. Apply projection
//
//     Add pressure gradient term back to u*:
//
//      u** = u* + dt * grad p / rho^nph
//
//     Solve Poisson equation for phi:
//
//     div( grad(phi) / rho^nph ) = div( u** )
//
//     Update pressure:
//
//     p = phi / dt
//
//     Update velocity, now divergence free
//
//     vel = u** - dt * grad p / rho^nph
//
// It is assumed that the ghost cels of the old data have been filled and
// the old and new data are the same in valid region.
//
void incflo::ApplyPredictor (bool incremental_projection)
{
    BL_PROFILE("amr-wind::incflo::ApplyPredictor")

    // We use the new time value for things computed on the "*" state
    Real new_time = m_time.new_time();

    if (m_verbose > 2) PrintMaxValues("before predictor step");

    if (m_use_godunov)
        amr_wind::io::print_mlmg_header("Godunov:");
    else
        amr_wind::io::print_mlmg_header("Predictor:");

    auto& icns_fields = icns().fields();
    auto& velocity_new = icns_fields.field;
    auto& velocity_old = velocity_new.state(amr_wind::FieldState::Old);
    auto& density_new = density();
    auto& density_old = density_new.state(amr_wind::FieldState::Old);
    auto& density_nph = density_new.state(amr_wind::FieldState::NPH);

    auto& velocity_forces = icns_fields.src_term;
    // only the old states are used in predictor
    auto& divtau = m_use_godunov
                       ? icns_fields.diff_term
                       : icns_fields.diff_term.state(amr_wind::FieldState::Old);

    // *************************************************************************************
    // Define the forcing terms to use in the Godunov prediction
    // *************************************************************************************
    if (m_use_godunov)
    {
        icns().compute_source_term(amr_wind::FieldState::Old);
        for (auto& seqn: scalar_eqns()) {
            seqn->compute_source_term(amr_wind::FieldState::Old);
        }
    }

    // *************************************************************************************
    // Compute viscosity / diffusive coefficients
    // *************************************************************************************
    m_sim.turbulence_model().update_turbulent_viscosity(amr_wind::FieldState::Old);
    icns().compute_mueff(amr_wind::FieldState::Old);
    for (auto& eqns: scalar_eqns())
        eqns->compute_mueff(amr_wind::FieldState::Old);

    // *************************************************************************************
    // Compute explicit viscous term
    // *************************************************************************************
    if (need_divtau()) {
        // Reuse existing buffer to avoid creating new multifabs
        amr_wind::field_ops::copy(velocity_new, velocity_old, 0, 0, velocity_new.num_comp(), 1);
        icns().compute_diffusion_term(amr_wind::FieldState::Old);
        if (m_use_godunov)
            amr_wind::field_ops::add(velocity_forces, divtau, 0, 0, AMREX_SPACEDIM, 0);
    }



    // *************************************************************************************
    // Compute explicit diffusive terms
    // *************************************************************************************
    if (need_divtau()) {
        for (auto& eqn: scalar_eqns()) {
            auto& field = eqn->fields().field;
            // Reuse existing buffer to avoid creating new multifabs
            amr_wind::field_ops::copy(field, field.state(amr_wind::FieldState::Old),
                                      0, 0, field.num_comp(), 1);

            eqn->compute_diffusion_term(amr_wind::FieldState::Old);

            if (m_use_godunov)
                amr_wind::field_ops::add(
                    eqn->fields().src_term,
                    eqn->fields().diff_term, 0, 0,
                    field.num_comp(), 0);
        }
    }

    if (m_use_godunov) {
       IntVect ng(nghost_force());
       icns().fields().src_term.fillpatch(m_time.current_time(), ng);

       for (auto& eqn: scalar_eqns()) {
           eqn->fields().src_term.fillpatch(m_time.current_time(), ng);
       }
    }

    // *************************************************************************************
    // if ( m_use_godunov) Compute the explicit advective terms
    //                     R_u^(n+1/2), R_s^(n+1/2) and R_t^(n+1/2)
    // if (!m_use_godunov) Compute the explicit advective terms
    //                     R_u^n      , R_s^n       and R_t^n
    // *************************************************************************************
    icns().compute_advection_term(amr_wind::FieldState::Old);
    for (auto& seqn: scalar_eqns()) {
        seqn->compute_advection_term(amr_wind::FieldState::Old);
    }

    // *************************************************************************************
    // Update density first
    // *************************************************************************************
    if (m_constant_density)
    {
        amr_wind::field_ops::copy(density_nph, density_old, 0, 0, 1, 1);
    }

    // Perform scalar update one at a time. This is to allow an updated density
    // at `n+1/2` to be computed before other scalars use it when computing
    // their source terms.
    for (auto& eqn: scalar_eqns()) {
        // Compute (recompute for Godunov) the scalar forcing terms
        eqn->compute_source_term(amr_wind::FieldState::NPH);

        // Update the scalar (if explicit), or the RHS for implicit/CN
        eqn->compute_predictor_rhs(m_diff_type);

        auto& field = eqn->fields().field;
        if (m_diff_type != DiffusionType::Explicit) {
            amrex::Real dt_diff = (m_diff_type == DiffusionType::Implicit)
                ? m_time.deltaT() : 0.5 * m_time.deltaT();

            // Solve diffusion eqn. and update of the scalar field
            eqn->solve(dt_diff);
        }

        // Update scalar at n+1/2
        amr_wind::field_ops::lincomb(
            field.state(amr_wind::FieldState::NPH),
            0.5, field.state(amr_wind::FieldState::Old), 0,
            0.5, field, 0, 0, field.num_comp(), 1);
    }

    // *************************************************************************************
    // Define (or if use_godunov, re-define) the forcing terms, without the viscous terms
    //    and using the half-time density
    // *************************************************************************************
    icns().compute_source_term(amr_wind::FieldState::New);

    // *************************************************************************************
    // Update the velocity
    // *************************************************************************************
    icns().compute_predictor_rhs(m_diff_type);

    // *************************************************************************************
    // Solve diffusion equation for u* but using eta_old at old time
    // *************************************************************************************
    if (m_diff_type == DiffusionType::Crank_Nicolson ||
        m_diff_type == DiffusionType::Implicit) {
        Real dt_diff = (m_diff_type == DiffusionType::Implicit)
                           ? m_time.deltaT()
                           : 0.5 * m_time.deltaT();
        icns().solve(dt_diff);
    }

    // ************************************************************************************
    //
    // Project velocity field, update pressure
    //
    // ************************************************************************************
    ApplyProjection(
        (density_nph).vec_const_ptrs(), new_time, m_time.deltaT(),
        incremental_projection);
}


//
// Apply corrector:
//
//  Output variables from the predictor are labelled _pred
//
//  1. Use u = vel_pred to compute
//
//      conv_u  = - u grad u
//      conv_r  = - u grad rho
//      conv_t  = - u grad trac
//      eta     = viscosity
//      divtau  = div( eta ( (grad u) + (grad u)^T ) ) / rho
//
//      conv_u  = 0.5 (conv_u + conv_u_pred)
//      conv_r  = 0.5 (conv_r + conv_r_pred)
//      conv_t  = 0.5 (conv_t + conv_t_pred)
//      if (m_diff_type == DiffusionType::Explicit)
//         divtau  = divtau at new_time using (*) state
//      else
//         divtau  = 0.0
//      eta     = eta at new_time
//
//     rhs = u + dt * ( conv + divtau )
//
//  2. Add explicit forcing term i.e. gravity + lagged pressure gradient
//
//      rhs += dt * ( g - grad(p + p0) / rho )
//
//      Note that in order to add the pressure gradient terms divided by rho,
//      we convert the velocity to momentum before adding and then convert them back.
//
//  3. A. If (m_diff_type == DiffusionType::Implicit)
//        solve implicit diffusion equation for u*
//
//     ( 1 - dt / rho * div ( eta grad ) ) u* = u^n + dt * conv_u
//                                                  + dt * ( g - grad(p + p0) / rho )
//
//     B. If (m_diff_type == DiffusionType::Crank-Nicolson)
//        solve semi-implicit diffusion equation for u*
//
//     ( 1 - (dt/2) / rho * div ( eta grad ) ) u* = u^n + dt * conv_u + (dt/2) / rho * div (eta_old grad) u^n
//                                                      + dt * ( g - grad(p + p0) / rho )
//
//  4. Apply projection
//
//     Add pressure gradient term back to u*:
//
//      u** = u* + dt * grad p / rho
//
//     Solve Poisson equation for phi:
//
//     div( grad(phi) / rho ) = div( u** )
//
//     Update pressure:
//
//     p = phi / dt
//
//     Update velocity, now divergence free
//
//     vel = u** - dt * grad p / rho
//
void incflo::ApplyCorrector()
{
    BL_PROFILE("amr-wind::incflo::ApplyCorrector")

    // We use the new time value for things computed on the "*" state
    Real new_time = m_time.new_time();

    if (m_verbose > 2) PrintMaxValues("before corrector step");

    amr_wind::io::print_mlmg_header("Corrector:");

    auto& density_new = density();
    auto& density_old = density_new.state(amr_wind::FieldState::Old);
    auto& density_nph = density_new.state(amr_wind::FieldState::NPH);

    // *************************************************************************************
    // Compute the explicit "new" advective terms R_u^(n+1,*), R_r^(n+1,*) and R_t^(n+1,*)
    // We only reach the corrector if !m_use_godunov which means we don't use the forces
    // in constructing the advection term
    // *************************************************************************************
    icns().compute_advection_term(amr_wind::FieldState::New);
    for (auto& seqn: scalar_eqns()) {
        seqn->compute_advection_term(amr_wind::FieldState::New);
    }

    // *************************************************************************************
    // Compute viscosity / diffusive coefficients
    // *************************************************************************************
    m_sim.turbulence_model().update_turbulent_viscosity(amr_wind::FieldState::New);
    icns().compute_mueff(amr_wind::FieldState::New);
    for (auto& eqns: scalar_eqns())
        eqns->compute_mueff(amr_wind::FieldState::New);

    // Here we create divtau of the (n+1,*) state that was computed in the predictor;
    //      we use this laps only if DiffusionType::Explicit
    if (m_diff_type == DiffusionType::Explicit) {
        icns().compute_diffusion_term(amr_wind::FieldState::New);

        for (auto& eqns: scalar_eqns()) {
            eqns->compute_diffusion_term(amr_wind::FieldState::New);
        }
    }

    // *************************************************************************************
    // Update density first
    // *************************************************************************************
    if (m_constant_density) {
        amr_wind::field_ops::copy(density_nph, density_old, 0, 0, 1, 1);
    }

    // Perform scalar update one at a time. This is to allow an updated density
    // at `n+1/2` to be computed before other scalars use it when computing
    // their source terms.
    for (auto& eqn: scalar_eqns()) {
        // Compute (recompute for Godunov) the scalar forcing terms
        // Note this is (rho * scalar) and not just scalar
        eqn->compute_source_term(amr_wind::FieldState::New);

        // Update (note that dtdt already has rho in it)
        // (rho trac)^new = (rho trac)^old + dt * (
        //                   div(rho trac u) + div (mu grad trac) + rho * f_t
        eqn->compute_corrector_rhs(m_diff_type);

        auto& field = eqn->fields().field;
        if (m_diff_type != DiffusionType::Explicit) {
            amrex::Real dt_diff = (m_diff_type == DiffusionType::Implicit)
                ? m_time.deltaT() : 0.5 * m_time.deltaT();

            // Solve diffusion eqn. and update of the scalar field
            eqn->solve(dt_diff);
        }

        // Update scalar at n+1/2
        amr_wind::field_ops::lincomb(
            field.state(amr_wind::FieldState::NPH),
            0.5, field.state(amr_wind::FieldState::Old), 0,
            0.5, field, 0, 0, field.num_comp(), 1);
    }

    // *************************************************************************************
    // Define the forcing terms to use in the final update (using half-time density)
    // *************************************************************************************
    icns().compute_source_term(amr_wind::FieldState::New);

    // *************************************************************************************
    // Update velocity
    // *************************************************************************************
    icns().compute_corrector_rhs(m_diff_type);

    // *************************************************************************************
    //
    // Solve diffusion equation for u* at t^{n+1} but using eta at predicted new time
    //
    // *************************************************************************************

    if (m_diff_type == DiffusionType::Crank_Nicolson ||
        m_diff_type == DiffusionType::Implicit) {
        Real dt_diff = (m_diff_type == DiffusionType::Implicit)
                           ? m_time.deltaT()
                           : 0.5 * m_time.deltaT();
        icns().solve(dt_diff);
    }

    // *************************************************************************************
    // Project velocity field, update pressure
    // *************************************************************************************
    bool incremental = false;
    ApplyProjection((density_nph).vec_const_ptrs(),new_time, m_time.deltaT(), incremental);

}
