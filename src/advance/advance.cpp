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

void incflo::Advance()
{
	BL_PROFILE("incflo::Advance");

    // Start timing current time step
    Real strt_step = ParallelDescriptor::second();

    if(verbose > 0) amrex::Print() << "\n ============   NEW TIME STEP   ============ \n";

    bool proj_2 = true;

    // Fill ghost nodes and reimpose boundary conditions
    for(int lev = 0; lev <= finest_level; lev++)
    {
        fill_mf_bc(lev, *ro[lev]);
        fill_mf_bc(lev, *eta[lev]);
    }
    incflo_set_velocity_bcs(t, 0);

    // Create temporary multifabs to hold the old-time conv and divtau
    // so we don't have to re-compute them in the corrector
    // TODO: Make this smarter -- actually shouldn't need separate implementatino of applyinf pred corr
    Vector<std::unique_ptr<MultiFab>> conv_old;
    Vector<std::unique_ptr<MultiFab>> divtau_old;
    conv_old.resize(finest_level + 1);
    divtau_old.resize(finest_level + 1);
    for(int lev = 0; lev <= finest_level; lev++)
    {
        conv_old[lev].reset(new MultiFab(grids[lev], dmap[lev], 3, 0, MFInfo(), *ebfactory[lev]));
        divtau_old[lev].reset(new MultiFab(grids[lev], dmap[lev], 3, 0, MFInfo(), *ebfactory[lev]));
    }

    // Compute time step size
    int initialisation = 0;
    incflo_compute_dt(initialisation);

    if(verbose > 0)
    {
        amrex::Print() << "\nStep " << nstep + 1 
                       << ": from old_time " << t 
                       << " to new time " << t + dt 
                       << " with dt = " << dt << ".\n" << std::endl;
    }

    // Backup velocity to old
    for(int lev = 0; lev <= finest_level; lev++)
        MultiFab::Copy(*vel_o[lev], *vel[lev], 0, 0, vel[lev]->nComp(), vel_o[lev]->nGrow());

    // Predictor step
    incflo_apply_predictor(conv_old, divtau_old, proj_2);

    // Corrector step
    incflo_apply_corrector(conv_old, divtau_old, proj_2);

    // Stop timing current time step
    Real end_step = ParallelDescriptor::second() - strt_step;
    ParallelDescriptor::ReduceRealMax(end_step, ParallelDescriptor::IOProcessorNumber());
    amrex::Print() << "Time per step " << end_step << std::endl;

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
void incflo::incflo_compute_dt(int initialisation)
{
	// DT is always computed even for fixed dt, so we can
	// issue a warning if fixed dt does not satisfy CFL condition.
	Real dt_new = dt;

	// Compute dt for this time step
	Real umax = -1.e20;
	Real vmax = -1.e20;
	Real wmax = -1.e20;
	Real romin = 1.e20;
	Real etamax = 0.0;

    // We only compute gp0max on the coarsest level because it is the same at all 
	Real gp0max[3];
    gp0max[0] = incflo_norm0(gp0, 0, 0);
    gp0max[1] = incflo_norm0(gp0, 0, 1);
    gp0max[2] = incflo_norm0(gp0, 0, 2);

    for(int lev = 0; lev <= finest_level; lev++)
    {
        umax = std::max(umax, incflo_norm0(vel, lev, 0));
        vmax = std::max(vmax, incflo_norm0(vel, lev, 1));
        wmax = std::max(wmax, incflo_norm0(vel, lev, 2));
        romin = std::min(romin, incflo_norm0(ro, lev, 0));
        // WARNING: This may cause trouble as we are not doing fully implicit solve!
        // TODO: revisit after testing fully 2/3 dimensional flows
        if(explicit_diffusion)
            etamax = std::max(etamax, incflo_norm0(eta, lev, 0));
    }

    const Real* dx = geom[finest_level].CellSize();
    Real idx = 1.0 / dx[0];
    Real idy = 1.0 / dx[1];
    Real idz = 1.0 / dx[2];

    // Convective term
    Real conv_cfl = umax * idx + vmax * idy + wmax * idz;

    // Viscous term
    Real diff_cfl = 2.0 * etamax / romin * (idx * idx + idy * idy + idz * idz);

    // Forcing term 
    Real forc_cfl = std::abs(gravity[0] - gp0max[0]) * idx 
                  + std::abs(gravity[1] - gp0max[1]) * idy  
                  + std::abs(gravity[2] - gp0max[2]) * idz;

    // Combined CFL conditioner 
    Real comb_cfl = conv_cfl + diff_cfl + sqrt(pow(conv_cfl + diff_cfl, 2) + 4.0 * forc_cfl);

    // Update dt
    dt_new = 2.0 * cfl / comb_cfl;

    // Reduce CFL for initial step
    if(initialisation)
        dt_new *= 0.1;

    // Protect against very small comb_cfl
    // This may happen, for example, when the initial velocity field
    // is zero for an inviscid flow with no external forcing
    if(comb_cfl <= 1.0e-18)
        dt_new = 0.5 * dt;

    // Don't let the timestep grow by more than 10% per step.
    if(dt > 0.0)
        dt_new = std::min(dt_new, 1.1*dt);

    // Don't overshoot the final time if not running to steady state
    if((!steady_state) & (stop_time > 0.0))
    {
        if(t + dt_new > stop_time)
            dt_new = stop_time - t;
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
// Compute predictor:
//
//  1. Compute
//
//     vel = vel_o + dt * R_u^n + dt * divtau*(1/ro)
//
//  2. Add explicit forcing term ( AKA gravity, lagged pressure gradient)
//
//     vel = vel + dt * ( g - grad(p+p0)/ro)
//
//  3. Add implicit forcing term 
//  
//     vel = vel / ( 1 + dt * f_gds/ro )
//
//  4. Solve for phi
//
//     div( grad(phi) / ro ) = div( vel / dt + grad(p)/ro )
//
//  5. Compute
//
//     vel = vel -  dt * grad(phi) / ro
//
//  6. Define
//
//     p = phi
//
void incflo::incflo_apply_predictor(Vector<std::unique_ptr<MultiFab>>& conv_old, 
                                    Vector<std::unique_ptr<MultiFab>>& divtau_old, 
                                    bool proj_2)
{
	BL_PROFILE("incflo::incflo_apply_predictor");

    // We use the new ime value for things computed on the "*" state
    Real new_time = t + dt; 

    // Compute the explicit advective term R_u^n
    incflo_compute_ugradu_predictor(conv_old, vel_o);

    update_derived_quantities();

    for(int lev = 0; lev <= finest_level; lev++)
    {
        // explicit_diffusion == true:  compute the full diffusive terms here
        // explicit_diffusion == false: compute only the off-diagonal terms here
        incflo_compute_divtau(lev, *divtau_old[lev], vel_o);

        // First add the convective term
        MultiFab::Saxpy(*vel[lev], dt, *conv_old[lev], 0, 0, 3, 0);

        // Add the diffusion terms (either all if explicit_diffusion == true or just
        // the off-diagonal terms if explicit_diffusion == false)
        MultiFab::Saxpy(*vel[lev], dt, *divtau_old[lev], 0, 0, 3, 0);

        // Add gravitational forces
        for(int dir = 0; dir < 3; dir++)
            (*vel[lev]).plus(dt * gravity[dir], dir, 1, 0);

        // Convert velocities to momenta
        for(int dir = 0; dir < 3; dir++)
            MultiFab::Multiply(*vel[lev], (*ro[lev]), 0, dir, 1, vel[lev]->nGrow());

        // Add (-dt grad p to momenta)
        MultiFab::Saxpy(*vel[lev], -dt, *gp[lev], 0, 0, 3, vel[lev]->nGrow());
        MultiFab::Saxpy(*vel[lev], -dt, *gp0[lev], 0, 0, 3, vel[lev]->nGrow());

        // Convert momenta back to velocities
        for(int dir = 0; dir < 3; dir++)
            MultiFab::Divide(*vel[lev], (*ro[lev]), 0, dir, 1, vel[lev]->nGrow());
    }

    // If doing implicit diffusion, solve here for u^*
    if(!explicit_diffusion) incflo_diffuse_velocity(new_time);

	// Project velocity field
	incflo_apply_projection(new_time, dt, proj_2);
    
	incflo_set_velocity_bcs(new_time, 0);

    if(verbose > 1)
    {
        amrex::Print() << "\nAfter predictor step:\n";
        for(int lev = 0; lev <= finest_level; lev++)
            incflo_print_max_vel(lev);

        incflo_compute_divu(t + dt);
        for(int lev = 0; lev <= finest_level; lev++)
            amrex::Print() << "max(abs(divu)) = " << incflo_norm0(divu, lev, 0) << "\n";
    }
}

//
// Compute corrector:
//
//  1. Compute
//
//     vel = vel_o + dt * (R_u^* + R_u^n) / 2 + dt * divtau*(1/ro)
//
//     where the starred variables are computed using "predictor-step"
//     variables.
//
//  2. Add explicit forcing term ( AKA gravity, lagged pressure gradient )
//
//     vel = vel + dt * ( g - grad(p+p0)/ro)
//
//  3. Add implicit forcing term 
//
//     vel = vel / ( 1 + dt * f_gds/ro )
//
//  4. Solve for phi
//
//     div( grad(phi) / ro ) = div(  vel / dt + grad(p)/ro )
//
//  5. Compute
//
//     vel = vel -  dt * grad(phi) / ro
//
//  6. Define
//
//     p = phi
//
void incflo::incflo_apply_corrector(Vector<std::unique_ptr<MultiFab>>& conv_old, 
                                    Vector<std::unique_ptr<MultiFab>>& divtau_old, 
                                    bool proj_2)
{
	BL_PROFILE("incflo::incflo_apply_corrector");

    // We use the new time value for things computed on the "*" state
    Real new_time = t + dt; 

	Vector<std::unique_ptr<MultiFab>> conv;
	Vector<std::unique_ptr<MultiFab>> divtau;   
    conv.resize(finest_level + 1);
    divtau.resize(finest_level + 1);
    for(int lev = 0; lev <= finest_level; lev++)
    {
        conv[lev].reset(new MultiFab(grids[lev], dmap[lev], 3, 0, MFInfo(), *ebfactory[lev]));
        divtau[lev].reset(new MultiFab(grids[lev], dmap[lev], 3, 0, MFInfo(), *ebfactory[lev]));
    }

    // Compute the explicit advective term R_u^*
    incflo_compute_ugradu_corrector(conv, vel);

    update_derived_quantities();

    for(int lev = 0; lev <= finest_level; lev++)
    {
        // If explicit_diffusion == true  then we compute the full diffusive terms here
        // If explicit_diffusion == false then we compute only the off-diagonal terms here
        incflo_compute_divtau(lev, *divtau[lev], vel);

        // Define u = u_o + dt/2 (R_u^* + R_u^n)
        MultiFab::LinComb(*vel[lev], 1.0, *vel_o[lev], 0, dt / 2.0, *conv[lev], 0, 0, 3, 0);
        MultiFab::Saxpy(*vel[lev], dt / 2.0, *conv_old[lev], 0, 0, 3, 0);

        // Add the diffusion terms (either all if explicit_diffusion == true or just
        // the
        //    off-diagonal terms if explicit_diffusion == false)
        MultiFab::Saxpy(*vel[lev], dt / 2.0, *divtau[lev], 0, 0, 3, 0);
        MultiFab::Saxpy(*vel[lev], dt / 2.0, *divtau_old[lev], 0, 0, 3, 0);

        // Add gravitational forces
        for(int dir = 0; dir < 3; dir++)
            (*vel[lev]).plus(dt * gravity[dir], dir, 1, 0);

        // Convert velocities to momenta
        for(int dir = 0; dir < 3; dir++)
            MultiFab::Multiply(*vel[lev], (*ro[lev]), 0, dir, 1, vel[lev]->nGrow());

        // Add (-dt grad p to momenta)
        MultiFab::Saxpy(*vel[lev], -dt, *gp[lev], 0, 0, 3, vel[lev]->nGrow());
        MultiFab::Saxpy(*vel[lev], -dt, *gp0[lev], 0, 0, 3, vel[lev]->nGrow());

        // Convert momenta back to velocities
        for(int dir = 0; dir < 3; dir++)
            MultiFab::Divide(*vel[lev], (*ro[lev]), 0, dir, 1, vel[lev]->nGrow());
    }

    // If doing implicit diffusion, solve here for u^*
    if(!explicit_diffusion) incflo_diffuse_velocity(new_time);

	// Apply projection
	incflo_apply_projection(new_time, dt, proj_2);
    
	// Project velocity field
	incflo_set_velocity_bcs(new_time, 0);

    if(verbose > 1)
    {
        amrex::Print() << "\nAfter corrector step:\n";
        for(int lev = 0; lev <= finest_level; lev++)
            incflo_print_max_vel(lev);

        incflo_compute_divu(t + dt);
        for(int lev = 0; lev <= finest_level; lev++)
            amrex::Print() << "max(abs(divu)) = " << incflo_norm0(divu, lev, 0) << "\n";
    }
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
int incflo::steady_state_reached()
{
    int condition1[finest_level + 1];
    int condition2[finest_level + 1];

    // Make sure velocity is up to date
    incflo_set_velocity_bcs(t, 0);

    for(int lev = 0; lev <= finest_level; lev++)
    {
        // Use temporaries to store the difference between current and previous solution
        MultiFab diff_vel(vel[lev]->boxArray(), vel[lev]->DistributionMap(), 3, 0);
        MultiFab::LinComb(diff_vel, 1.0, *vel[lev], 0, -1.0, *vel_o[lev], 0, 0, 3, 0);

        Real max_change = 0.0;
        Real max_relchange = 0.0;
        // Loop over components, only need to check the largest one
        for(int i = 0; i < 3; i++)
        {
            // max(abs(u^{n+1}-u^n))
            max_change = std::max(max_change, incflo_norm0(diff_vel, lev, i));
            
            // sum(abs(u^{n+1}-u^n)) / sum(abs(u^n))
            Real norm1_diff = incflo_norm1(diff_vel, lev, i);
            Real norm1_old = incflo_norm1(vel_o, lev, i);
            Real relchange = norm1_old > 1.0e-15 ? norm1_diff / norm1_old : 0.0;
            max_relchange = std::max(max_relchange, relchange);
        }

        condition1[lev] = (max_change < steady_state_tol * dt);
        condition2[lev] = (max_relchange < steady_state_tol);

        // Print out info on steady state checks
        if(verbose > 0)
        {
            amrex::Print() << "\nSteady state check:\n";
            amrex::Print() << "||u-uo||/||uo|| = " << max_relchange 
                           << ", du/dt  = " << max_change/dt << std::endl;
        }
    }

    int reached = 1;
    for(int lev = 0; lev <= finest_level; lev++)
    {
        reached = reached && (condition1[lev] || condition2[lev]);
    }

	// Always return negative to first access. This way
	// initial zero velocity field do not test for false positive
	if(nstep < 2)
		return 0;
	else
		return reached;
}
