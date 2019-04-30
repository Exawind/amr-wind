#include <AMReX_Array.H>
#include <AMReX_BC_TYPES.H>
#include <AMReX_BLassert.H>
#include <AMReX_Box.H>
#include <AMReX_EBMultiFabUtil.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_VisMF.H>

#include <incflo.H>
#include <mac_F.H>
#include <projection_F.H>
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
    FillScalarBC();
    FillVelocityBC(cur_time, 0);

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
        MultiFab::Copy(*vel_o[lev], *vel[lev], 0, 0, vel[lev]->nComp(), vel_o[lev]->nGrow());
    }

    ApplyPredictor();

    ApplyCorrector();

    if(incflo_verbose > 1)
    {
        amrex::Print() << "End of time step: " << std::endl;
        PrintMaxValues(cur_time + dt);
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
// Compute new dt by using the formula derived in
// "A Boundary Condition Capturing Method for Multiphase Incompressible Flow"
// by Kang et al. (JCP).
//
//  dt/2 * ( C+V + sqrt( (C+V)**2 + 4Fx/dx + 4Fy/dy + 4Fz/dz )
//
// where
//
// C = max(|U|)/dx + max(|V|)/dy + max(|W|)/dz    --> Convection
//
// V = 2 * max(eta/ro) * (1/dx^2 + 1/dy^2 +1/dz^2) --> Diffusion
//
// Fx, Fy, Fz = net acceleration due to external forces
//
// WARNING: We use a slightly modified version of C in the implementation below
//
void incflo::ComputeDt(int initialisation)
{
	BL_PROFILE("incflo::ComputeDt");

	// Compute dt for this time step
	Real umax = 0.0;
	Real vmax = 0.0;
	Real wmax = 0.0;
	Real romin = 1.e20;
	Real etamax = 0.0;

    for(int lev = 0; lev <= finest_level; lev++)
    {
        // The functions take the min/max over uncovered cells 
        umax   = amrex::max(umax,   Norm(vel, lev, 0, 0));
        vmax   = amrex::max(vmax,   Norm(vel, lev, 1, 0));
        wmax   = amrex::max(wmax,   Norm(vel, lev, 2, 0));
        romin  = amrex::min(romin,  Norm( ro, lev, 0, 0));
        etamax = amrex::max(etamax, Norm(eta, lev, 0, 0));
    }

    const Real* dx = geom[finest_level].CellSize();
    Real idx = 1.0 / dx[0];
    Real idy = 1.0 / dx[1];
    Real idz = 1.0 / dx[2];

    // Convective term
    Real conv_cfl = std::max(std::max(umax * idx, vmax * idy), wmax * idz);

    // Viscous term
    Real diff_cfl = 2.0 * etamax / romin * (idx * idx + idy * idy + idz * idz);

    // Forcing term
    Real forc_cfl = std::abs(gravity[0] - std::abs(gp0[0])) * idx
                  + std::abs(gravity[1] - std::abs(gp0[1])) * idy
                  + std::abs(gravity[2] - std::abs(gp0[2])) * idz;

    // Combined CFL conditioner
    Real comb_cfl = conv_cfl + diff_cfl + sqrt(pow(conv_cfl + diff_cfl, 2) + 4.0 * forc_cfl);

    // Update dt
    Real dt_new = 2.0 * cfl / comb_cfl;

    // Reduce CFL for initial step
    if(initialisation)
    {
        dt_new *= 0.1;
    }

    // Protect against very small comb_cfl
    // This may happen, for example, when the initial velocity field
    // is zero for an inviscid flow with no external forcing
    Real eps = std::numeric_limits<Real>::epsilon();
    if(comb_cfl <= eps)
    {
        dt_new = 0.5 * dt;
    }

    // Don't let the timestep grow by more than 10% per step.
    if(dt > 0.0 && last_plt != nstep)
    {
        dt_new = amrex::min(dt_new, 1.1 * dt);
    }
    
    // Don't overshoot specified plot times
    if(plot_per > 0.0 && 
            (trunc((cur_time + dt_new + eps) / plot_per) > trunc((cur_time + eps) / plot_per)))
    {
        dt_new = trunc((cur_time + dt_new) / plot_per) * plot_per - cur_time;
    }

    // Don't overshoot the final time if not running to steady state
    if(!steady_state && stop_time > 0.0)
    {
        if(cur_time + dt_new > stop_time)
        {
            dt_new = stop_time - cur_time;
        }
    }

    // Make sure the timestep is not set to zero after a plot_per stop
    if(dt_new < eps)
    {
        dt_new = 0.5 * dt;
    }

    // If using fixed time step, check CFL condition and give warning if not satisfied
	if(fixed_dt > 0.0)
	{
		if(dt_new < fixed_dt)
		{
			amrex::Print() << "WARNING: fixed_dt does not satisfy CFL condition: \n"
						   << "max dt by CFL     : " << dt_new << "\n"
						   << "fixed dt specified: " << fixed_dt << std::endl;
		}
		dt = fixed_dt;
	}
	else
	{
		dt = dt_new;
	}
}

//
// Apply predictor:
//
//  1. Use u = vel_old to compute
//
//      conv    = - u grad u
//      eta     = eta( ||strainrate|| ) 
//      divtau  = div( eta (grad u)^T ) / rho
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
//     div( grad(phi) / ro ) = div( u** )
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

    if(incflo_verbose > 2)
    {
        amrex::Print() << "Before predictor step:" << std::endl;
        PrintMaxValues(new_time);
    }

    // Compute the explicit advective term: conv = - u dot grad(u)
    ComputeUGradU(conv_old, vel_o, cur_time);

    // Update the derived quantities, notably strain-rate tensor and viscosity
    UpdateDerivedQuantities();

    for(int lev = 0; lev <= finest_level; lev++)
    {
        // Save this value of eta as eta_old for use in the corrector as well
        MultiFab::Copy(*eta_old[lev], *eta[lev], 0, 0, eta[lev]->nComp(), eta_old[lev]->nGrow());

        // compute only the off-diagonal terms here
        ComputeDivTau(lev, *divtau_old[lev], vel_o);

        // First add the convective term
        MultiFab::Saxpy(*vel[lev], dt, *conv_old[lev], 0, 0, 3, 0);

        // Add the viscous terms         
        MultiFab::Saxpy(*vel[lev], dt, *divtau_old[lev], 0, 0, 3, 0);

        // Add gravitational forces
        for(int dir = 0; dir < 3; dir++)
        {
            (*vel[lev]).plus(dt * gravity[dir], dir, 1, 0);
        }

        // Convert velocities to momenta
        for(int dir = 0; dir < 3; dir++)
        {
            MultiFab::Multiply(*vel[lev], (*ro[lev]), 0, dir, 1, vel[lev]->nGrow());
        }

        // Add (-dt grad p to momenta)
        MultiFab::Saxpy(*vel[lev], -dt, *gp[lev], 0, 0, 3, vel[lev]->nGrow());
        for(int dir = 0; dir < 3; dir++)
        {
            (*vel[lev]).plus(-dt * gp0[dir], dir, 1, 0);
        }

        // Convert momenta back to velocities
        for(int dir = 0; dir < 3; dir++)
        {
            MultiFab::Divide(*vel[lev], (*ro[lev]), 0, dir, 1, vel[lev]->nGrow());
        }
    }
    FillVelocityBC(new_time, 0);

    // Solve implicit diffusion equation for u*
    diffusion_equation->solve(vel, ro, eta, dt);

	// Project velocity field, update pressure
	ApplyProjection(new_time, dt);

    // Fill velocity BCs again
	FillVelocityBC(new_time, 0);
}

//
// Apply corrector:
//
//  Output variables from the predictor are labelled _pred 
//
//  1. Use u = vel_pred to compute
//
//      conv    = - u grad u
//      eta     = eta( ||strainrate|| ) 
//      divtau  = div( eta (grad u)^T ) / rho
//
//      conv    = 0.5 (conv + conv_pred)
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
//     div( grad(phi) / ro ) = div( u** )
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

    // Compute the explicit advective term: conv = - u dot grad(u)
    ComputeUGradU(conv, vel, new_time);

    // Update the derived quantities, notably strain-rate tensor and viscosity
    UpdateDerivedQuantities();

    for(int lev = 0; lev <= finest_level; lev++)
    {
        // compute only the off-diagonal terms here
        ComputeDivTau(lev, *divtau[lev], vel);

        // First add the convective terms
        MultiFab::LinComb(*vel[lev], 1.0, *vel_o[lev], 0, dt / 2.0, *conv[lev], 0, 0, 3, 0);
        MultiFab::Saxpy(*vel[lev], dt / 2.0, *conv_old[lev], 0, 0, 3, 0);

        // Add the viscous terms         
        MultiFab::Saxpy(*vel[lev], dt / 2.0, *divtau[lev], 0, 0, 3, 0);
        MultiFab::Saxpy(*vel[lev], dt / 2.0, *divtau_old[lev], 0, 0, 3, 0);

        // Add gravitational forces
        for(int dir = 0; dir < 3; dir++)
        {
            (*vel[lev]).plus(dt * gravity[dir], dir, 1, 0);
        }

        // Convert velocities to momenta
        for(int dir = 0; dir < 3; dir++)
        {
            MultiFab::Multiply(*vel[lev], (*ro[lev]), 0, dir, 1, vel[lev]->nGrow());
        }

        // Add (-dt grad p to momenta)
        MultiFab::Saxpy(*vel[lev], -dt, *gp[lev], 0, 0, 3, vel[lev]->nGrow());
        for(int dir = 0; dir < 3; dir++)
        {
            (*vel[lev]).plus(-dt * gp0[dir], dir, 1, 0);
        }

        // Convert momenta back to velocities
        for(int dir = 0; dir < 3; dir++)
        {
            MultiFab::Divide(*vel[lev], (*ro[lev]), 0, dir, 1, vel[lev]->nGrow());
        }

        // Take eta as the average of the predictor and corrector values
        MultiFab::LinComb(*eta[lev], 0.5, *eta_old[lev], 0, 0.5, *eta[lev], 0, 0, 1, 0);
    }
    FillVelocityBC(new_time, 0);

    // Solve implicit diffusion equation for u*
    diffusion_equation->solve(vel, ro, eta, dt);

	// Project velocity field, update pressure
	ApplyProjection(new_time, dt);

    // Fill velocity BCs again
	FillVelocityBC(new_time, 0);
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
    FillVelocityBC(cur_time, 0);

    // Use temporaries to store the difference between current and previous solution
	Vector<std::unique_ptr<MultiFab>> diff_vel;
    diff_vel.resize(finest_level + 1);
    for(int lev = 0; lev <= finest_level; lev++)
    {
        diff_vel[lev].reset(new MultiFab(grids[lev], dmap[lev], 3, 0, MFInfo(), *ebfactory[lev]));
        MultiFab::LinComb(*diff_vel[lev], 1.0, *vel[lev], 0, -1.0, *vel_o[lev], 0, 0, 3, 0);

        Real max_change = 0.0;
        Real max_relchange = 0.0;
        // Loop over components, only need to check the largest one
        for(int i = 0; i < 3; i++)
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
    }
    else
    {
        return reached;
    }
}
