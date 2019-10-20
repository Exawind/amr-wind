#include <AMReX_Array.H>
#include <AMReX_BC_TYPES.H>
#include <AMReX_BLassert.H>
#include <AMReX_Box.H>
#include <AMReX_MultiFab.H>
#include <AMReX_VisMF.H>

#include <incflo.H>
#include <incflo_proj_F.H>
#include <setup_F.H>

#include <limits>

void incflo::Advance()
{
    BL_PROFILE("incflo::Advance");

    // Start timing current time step
    Real strt_step = ParallelDescriptor::second();

    if(incflo_verbose > 0)
    {
        amrex::Print() << "\n ============   NEW TIME STEP   ============ \n";
    }

    // Fill ghost nodes and reimpose boundary conditions
    if (!constant_density)
       incflo_set_density_bcs(cur_time, density);
    if (advect_tracer)
       incflo_set_tracer_bcs(cur_time, tracer);
    incflo_set_velocity_bcs(cur_time, vel, 0);

    // Compute time step size
    int initialisation = 0;
    ComputeDt(initialisation);

    // Set new and old time to correctly use in fillpatching
    for(int lev = 0; lev <= finest_level; lev++)
    {
        t_old[lev] = cur_time; 
        t_new[lev] = cur_time + dt; 
    }

    if(incflo_verbose > 0)
    {
        amrex::Print() << "\nStep " << nstep + 1
                       << ": from old_time " << cur_time
                       << " to new time " << cur_time + dt
                       << " with dt = " << dt << ".\n" << std::endl;
    }

    // Backup velocity to old
    for(int lev = 0; lev <= finest_level; lev++)
    {
        MultiFab::Copy(    *vel_o[lev],     *vel[lev], 0, 0,     vel[lev]->nComp(),     vel_o[lev]->nGrow());
        MultiFab::Copy(*density_o[lev], *density[lev], 0, 0, density[lev]->nComp(), density_o[lev]->nGrow());
        MultiFab::Copy(* tracer_o[lev],  *tracer[lev], 0, 0,  tracer[lev]->nComp(),  tracer_o[lev]->nGrow());
    }

    ApplyPredictor();

    ApplyCorrector();

    if(incflo_verbose > 1)
    {
        amrex::Print() << "End of time step: " << std::endl;
        PrintMaxValues(cur_time + dt);
        if(probtype%10 == 3 or probtype == 5)
        {
            ComputeDrag();
            amrex::Print() << "Drag force = " << (*drag[0]).sum(0, false) << std::endl; 
        }
    }

    // Stop timing current time step
    Real end_step = ParallelDescriptor::second() - strt_step;
    ParallelDescriptor::ReduceRealMax(end_step, ParallelDescriptor::IOProcessorNumber());
    if(incflo_verbose > 0)
    {
        amrex::Print() << "Time per step " << end_step << std::endl;
    }

	BL_PROFILE_REGION_STOP("incflo::Advance");
}

//
// Apply predictor:
//
//  1. Use u = vel_old to compute
//
//      conv_u  = - u grad u
//      conv_r  = - div( u rho  )
//      conv_t  = - div( u trac )
//      eta     = visosity
//      divtau  = div( eta ( (grad u) + (grad u)^T ) ) / rho
//
//      rhs = u + dt * ( conv + divtau )
//
//  2. Add explicit forcing term i.e. gravity + lagged pressure gradient
//
//      rhs += dt * ( g - grad(p + p0) / rho )
//
//      Note that in order to add the pressure gradient terms divided by rho, 
//      we convert the velocity to momentum before adding and then convert them back. 
//
//  3. Solve implicit diffusion equation for u* 
//
//     ( 1 - dt / rho * div ( eta grad ) ) u* = rhs
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
void incflo::ApplyPredictor()
{
    BL_PROFILE("incflo::ApplyPredictor");

    // We use the new ime value for things computed on the "*" state
    Real new_time = cur_time + dt;

    incflo_set_velocity_bcs(cur_time, vel_o, 0);

    if(incflo_verbose > 2)
    {
        amrex::Print() << "Before predictor step:" << std::endl;
        PrintMaxValues(new_time);
    }

    // Compute the explicit advective terms R_u^n and R_s^n
    incflo_compute_convective_term( conv_u_old, conv_r_old, conv_t_old, vel_o, density_o, tracer_o, cur_time );

    // This fills the eta array (if non-Newtonian, then using strain-rate of velocity at time "cur_time")
    // We shouldn't need to do this again because the cur_time velocity is vel_o which hasn't changed
    // ComputeViscosity();

    for(int lev = 0; lev <= finest_level; lev++)
    {
        // Save this value of eta as eta_old for use in the corrector as well
        MultiFab::Copy(*eta_old[lev], *eta[lev], 0, 0, eta[lev]->nComp(), eta_old[lev]->nGrow());
    }

    // Compute explicit diffusion if used
    if (explicit_diffusion)
       ComputeDivTau(divtau_old, vel_o, density_o, eta_old);
    else
       for (int lev = 0; lev <= finest_level; lev++)
          divtau_old[lev]->setVal(0.);

    for(int lev = 0; lev <= finest_level; lev++)
    {
        // First add the convective term to the velocity
        MultiFab::Saxpy(*vel[lev], dt, *conv_u_old[lev], 0, 0, AMREX_SPACEDIM, 0);

        //   Now add the convective term to the density and tracer
        if (!constant_density)
           MultiFab::Saxpy(*density[lev], dt, *conv_r_old[lev], 0, 0,      1, 0);
        if (advect_tracer)
           MultiFab::Saxpy( *tracer[lev], dt, *conv_t_old[lev],  0, 0, ntrac, 0);

        // Add the viscous terms         
        if (explicit_diffusion)
            MultiFab::Saxpy(*vel[lev], dt, *divtau_old[lev], 0, 0, AMREX_SPACEDIM, 0);

        // Add gravitational forces
        if (use_boussinesq)
        {
           // This uses a Boussinesq approximation where the buoyancy depends on first tracer
           //      rather than density
           for(int dir = 0; dir < AMREX_SPACEDIM; dir++)
              MultiFab::Saxpy(*vel[lev], dt*gravity[dir], *tracer[lev], 0, dir, 1, 0);
        } else {
           for(int dir = 0; dir < AMREX_SPACEDIM; dir++)
               (*vel[lev]).plus(dt * gravity[dir], dir, 1, 0);
        }

        // Convert velocities to momenta
        for(int dir = 0; dir < AMREX_SPACEDIM; dir++)
        {
            MultiFab::Multiply(*vel[lev], (*density[lev]), 0, dir, 1, vel[lev]->nGrow());
        }

        // Add (-dt grad p to momenta)
        MultiFab::Saxpy(*vel[lev], -dt, *gp[lev], 0, 0, AMREX_SPACEDIM, vel[lev]->nGrow());
        for(int dir = 0; dir < AMREX_SPACEDIM; dir++)
        {
            (*vel[lev]).plus(-dt * gp0[dir], dir, 1, 0);
        }

        // Convert momenta back to velocities
        for(int dir = 0; dir < AMREX_SPACEDIM; dir++)
        {
            MultiFab::Divide(*vel[lev], (*density[lev]), 0, dir, 1, vel[lev]->nGrow());
        }
    }
    if (!constant_density)
       incflo_set_density_bcs(new_time, density);
    if (advect_tracer)
       incflo_set_tracer_bcs(new_time, tracer);
    incflo_set_velocity_bcs(new_time, vel, 0);

    // Solve implicit diffusion equation for u*
    if (!explicit_diffusion)
       diffusion_equation->solve(vel, density, eta, dt);

    // Project velocity field, update pressure
    ApplyProjection(new_time, dt);

    // Fill velocity BCs again
    incflo_set_velocity_bcs(new_time, vel, 0);
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
//      divtau  = 0.5 (divtau + divtau_pred)
//      eta     = 0.5 (eta + eta_pred)
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
//  3. Solve implicit diffusion equation for u* 
//
//     ( 1 - dt / rho * div ( eta grad ) ) u* = rhs
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
	BL_PROFILE("incflo::ApplyCorrector");

    // We use the new time value for things computed on the "*" state
    Real new_time = cur_time + dt;

    if(incflo_verbose > 2)
    {
        amrex::Print() << "Before corrector step:" << std::endl;
        PrintMaxValues(new_time);
    }

    // Compute the explicit advective terms R_u^* and R_s^*
    incflo_compute_convective_term( conv_u, conv_r, conv_t, vel, density, tracer, new_time );

    // This fills the eta array (if non-Newtonian, then using strain-rate of velocity at time "cur_time")
    ComputeViscosity();

    // Compute explicit diffusion if used -- note that even though we call this "explicit",
    //   the diffusion term does end up being time-centered so formally second-order
    if (explicit_diffusion)
       ComputeDivTau(divtau, vel, density, eta);
    else
       for (int lev = 0; lev <= finest_level; lev++)
          divtau[lev]->setVal(0.);

    for(int lev = 0; lev <= finest_level; lev++)
    {
        // First add the convective terms to velocity
        MultiFab::LinComb(*vel[lev], 1.0, *vel_o[lev], 0, dt / 2.0, *conv_u[lev], 0, 0, AMREX_SPACEDIM, 0);
        MultiFab::Saxpy(*vel[lev], dt / 2.0, *conv_u_old[lev], 0, 0, AMREX_SPACEDIM, 0);

        //   Now add the convective terms to density
        if (!constant_density)
        {
           MultiFab::LinComb(*density[lev], 1.0, *density_o[lev], 0, dt / 2.0, *conv_r[lev], 0, 0, 1, 0);
           MultiFab::Saxpy(*density[lev], dt / 2.0, *conv_r_old[lev], 0, 0, 1, 0);
        }

        //   Now add the convective terms to tracer
        if (advect_tracer)
        {
           MultiFab::LinComb(*tracer[lev], 1.0, *tracer_o[lev], 0, dt / 2.0, *conv_t[lev]    , 0, 0, tracer[lev]->nComp(), 0);
           MultiFab::Saxpy(  *tracer[lev], dt / 2.0                        , *conv_t_old[lev], 0, 0, tracer[lev]->nComp(), 0);
        }

        // Add the viscous terms         
        if (explicit_diffusion)
        {
            MultiFab::Saxpy(*vel[lev], dt / 2.0, *divtau[lev]    , 0, 0, AMREX_SPACEDIM, 0);
            MultiFab::Saxpy(*vel[lev], dt / 2.0, *divtau_old[lev], 0, 0, AMREX_SPACEDIM, 0);
        }

        // Add gravitational forces
        if (use_boussinesq)
        {
           // This uses a Boussinesq approximation where the buoyancy depends on tracer
           //      rather than density
           for(int dir = 0; dir < AMREX_SPACEDIM; dir++)
              MultiFab::Saxpy(*vel[lev], dt*gravity[dir], *tracer[lev], 0, dir, 1, 0);
        } else {
           for(int dir = 0; dir < AMREX_SPACEDIM; dir++)
            (*vel[lev]).plus(dt * gravity[dir], dir, 1, 0);
        }

        // Convert velocities to momenta
        for(int dir = 0; dir < AMREX_SPACEDIM; dir++)
        {
            MultiFab::Multiply(*vel[lev], (*density[lev]), 0, dir, 1, vel[lev]->nGrow());
        }

        // Add (-dt grad p to momenta)
        MultiFab::Saxpy(*vel[lev], -dt, *gp[lev], 0, 0, AMREX_SPACEDIM, vel[lev]->nGrow());
        for(int dir = 0; dir < AMREX_SPACEDIM; dir++)
        {
            (*vel[lev]).plus(-dt * gp0[dir], dir, 1, 0);
        }

        // Convert momenta back to velocities
        for(int dir = 0; dir < AMREX_SPACEDIM; dir++)
        {
            MultiFab::Divide(*vel[lev], (*density[lev]), 0, dir, 1, vel[lev]->nGrow());
        }

        // Take eta as the average of the predictor and corrector values
        MultiFab::LinComb(*eta[lev], 0.5, *eta_old[lev], 0, 0.5, *eta[lev], 0, 0, 1, 0);
    }

    if (!constant_density)
       incflo_set_density_bcs(new_time, density);
    if (advect_tracer)
       incflo_set_tracer_bcs(new_time, tracer);
    incflo_set_velocity_bcs(new_time, vel, 0);

    // Solve implicit diffusion equation for u*
    if (!explicit_diffusion)
       diffusion_equation->solve(vel, density, eta, dt);

    // Project velocity field, update pressure
    ApplyProjection(new_time, dt);

    // Fill velocity BCs again
    incflo_set_velocity_bcs(new_time, vel, 0);
}

//
// Check if steady state has been reached by verifying that
//
//      max(abs( u^(n+1) - u^(n) )) / dt < tol
//      max(abs( v^(n+1) - v^(n) )) / dt < tol
//      max(abs( w^(n+1) - w^(n) )) / dt < tol
//
//      OR
//
//      sum(abs( u^(n+1) - u^(n) )) / sum(abs( u^(n) )) < tol
//      sum(abs( v^(n+1) - v^(n) )) / sum(abs( v^(n) )) < tol
//      sum(abs( w^(n+1) - w^(n) )) / sum(abs( w^(n) )) < tol
//
bool incflo::SteadyStateReached()
{
    BL_PROFILE("incflo::SteadyStateReached()");

    int condition1[finest_level + 1];
    int condition2[finest_level + 1];

    // Make sure velocity is up to date
    incflo_set_velocity_bcs(cur_time, vel, 0);

    // Use temporaries to store the difference between current and previous solution
	Vector<std::unique_ptr<MultiFab>> diff_vel;
    diff_vel.resize(finest_level + 1);
    for(int lev = 0; lev <= finest_level; lev++)
    {
        diff_vel[lev].reset(new MultiFab(grids[lev], dmap[lev], AMREX_SPACEDIM, 0, MFInfo(), *ebfactory[lev]));
        MultiFab::LinComb(*diff_vel[lev], 1.0, *vel[lev], 0, -1.0, *vel_o[lev], 0, 0, AMREX_SPACEDIM, 0);

        Real max_change = 0.0;
        Real max_relchange = 0.0;
        // Loop over components, only need to check the largest one
        for(int i = 0; i < AMREX_SPACEDIM; i++)
        {
            // max(abs(u^{n+1}-u^n))
            max_change = amrex::max(max_change, Norm(diff_vel, lev, i, 0));

            // sum(abs(u^{n+1}-u^n)) / sum(abs(u^n))
            // TODO: this gives zero often, check for bug
            Real norm1_diff = Norm(diff_vel, lev, i, 1);
            Real norm1_old = Norm(vel_o, lev, i, 1);
            Real relchange = norm1_old > 1.0e-15 ? norm1_diff / norm1_old : 0.0;
            max_relchange = amrex::max(max_relchange, relchange);
        }

        condition1[lev] = (max_change < steady_state_tol * dt);
        condition2[lev] = (max_relchange < steady_state_tol);

        // Print out info on steady state checks
        if(incflo_verbose > 0)
        {
            amrex::Print() << "\nSteady state check level " << lev << std::endl; 
            amrex::Print() << "||u-uo||/||uo|| = " << max_relchange
                           << ", du/dt  = " << max_change/dt << std::endl;
        }
    }

    bool reached = true;
    for(int lev = 0; lev <= finest_level; lev++)
    {
        reached = reached && (condition1[lev] || condition2[lev]);
    }

    // Always return negative to first access. This way
    // initial zero velocity field do not test for false positive
    if(nstep < 2)
    {
        return false;
    } else {
        return reached;
    }
}
