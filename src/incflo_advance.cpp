#include <incflo.H>

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
    bool explicit_diffusion = (m_diff_type == DiffusionType::Explicit);
    ComputeDt(initialisation, explicit_diffusion);

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

    if (!use_godunov)
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
//      eta_old     = visosity at cur_time
//      if (m_diff_type == DiffusionType::Explicit)
//         divtau _old = div( eta ( (grad u) + (grad u)^T ) ) / rho
//         rhs = u + dt * ( conv + divtau_old )
//      else
//         divtau_old  = 0.0 
//         rhs = u + dt * conv 
//
//      eta     = eta at new_time 
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
//     ( 1 - (dt/2) / rho * div ( eta_old grad ) ) u* = u^n + dt * conv_u + (dt/2) / rho * div (eta_old grad) u^n
//                                                          + dt * ( g - grad(p + p0) / rho )
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

    // We use the new time value for things computed on the "*" state
    Real new_time = cur_time + dt;

    incflo_set_velocity_bcs(cur_time, vel_o, 0);

    if(incflo_verbose > 2)
    {
        amrex::Print() << "Before predictor step:" << std::endl;
        PrintMaxValues(new_time);
    }

    // if ( use_godunov) Compute the explicit advective terms R_u^(n+1/2), R_s^(n+1/2) and R_t^(n+1/2)
    // if (!use_godunov) Compute the explicit advective terms R_u^n      , R_s^n       and R_t^n
    incflo_compute_convective_term( conv_u_old, conv_r_old, conv_t_old, vel_o, density_o, tracer_o, cur_time );

    // This fills the eta_old array (if non-Newtonian, then using strain-rate of velocity at time "cur_time")
    ComputeViscosity(eta_old, cur_time);
    
    // Adds additional viscosity if LES is on
    add_eddy_viscosity(eta, eta_tracer, cur_time);
    
    // Compute explicit diffusion if used
    if (m_diff_type == DiffusionType::Explicit ||
        m_diff_type == DiffusionType::Crank_Nicolson)
    {
       int extrap_dir_bcs = 0;
       incflo_set_velocity_bcs (cur_time, vel_o, extrap_dir_bcs);
       diffusion_op->ComputeDivTau(divtau_old,    vel_o, density_o, eta_old);

       diffusion_op->ComputeLapS  (  laps_old, tracer_o, density_o, mu_s);
    } else {
       for (int lev = 0; lev <= finest_level; lev++)
       {
          divtau_old[lev]->setVal(0.);
            laps_old[lev]->setVal(0.);
       }
    }

    for (int lev = 0; lev <= finest_level; lev++)
    {
        // First add the convective term to the velocity
        MultiFab::Saxpy(*vel[lev], dt, *conv_u_old[lev], 0, 0, AMREX_SPACEDIM, 0);

        //   Now add the convective term to the density and tracer
        if (!constant_density)
            MultiFab::Saxpy(*density[lev], dt, *conv_r_old[lev], 0, 0,      1, 0);

        if (advect_tracer)
        {
            MultiFab::Saxpy( *tracer[lev], dt, *conv_t_old[lev],  0, 0, ntrac, 0);
            for (int i = 0; i < ntrac; i++)
            {
                if (m_diff_type == DiffusionType::Explicit) 
                    MultiFab::Saxpy(*tracer[lev],     dt, *laps_old[lev], i, i, 1, 0);
                else if (m_diff_type == DiffusionType::Crank_Nicolson)
                    MultiFab::Saxpy(*tracer[lev], 0.5*dt, *laps_old[lev], i, i, 1, 0);
            }
        }

        // Add the viscous terms         
        if (m_diff_type == DiffusionType::Explicit) 
            MultiFab::Saxpy(*vel[lev], dt, *divtau_old[lev], 0, 0, AMREX_SPACEDIM, 0);
        else if (m_diff_type == DiffusionType::Crank_Nicolson)
            MultiFab::Saxpy(*vel[lev], 0.5*dt, *divtau_old[lev], 0, 0, AMREX_SPACEDIM, 0);

        // Add gravitational forces
        if (!use_boussinesq)
        {
           for(int dir = 0; dir < AMREX_SPACEDIM; dir++)
               (*vel[lev]).plus(dt * gravity[dir], dir, 1, 0);
        }

        add_abl_source_terms(*vel[lev],*tracer[lev],vel[lev]->nGrow());
        
        // Convert velocities to momenta
        for(int dir = 0; dir < AMREX_SPACEDIM; dir++)
            MultiFab::Multiply(*vel[lev], (*density[lev]), 0, dir, 1, vel[lev]->nGrow());

        // Add (-dt grad p to momenta)
        MultiFab::Saxpy(*vel[lev], -dt, *gp[lev], 0, 0, AMREX_SPACEDIM, vel[lev]->nGrow());

        // Add (-dt grad p0 to momenta if not Boussinesq)
        if (!use_boussinesq)
            for(int dir = 0; dir < AMREX_SPACEDIM; dir++)
                (*vel[lev]).plus(-dt * gp0[dir], dir, 1, 0);

        // Convert momenta back to velocities
        for(int dir = 0; dir < AMREX_SPACEDIM; dir++)
            MultiFab::Divide(*vel[lev], (*density[lev]), 0, dir, 1, vel[lev]->nGrow());
    }

    if (!constant_density)
        incflo_set_density_bcs(new_time, density);
    if (advect_tracer)
        incflo_set_tracer_bcs(new_time, tracer);
    incflo_set_velocity_bcs(new_time, vel, 0);

    // Solve diffusion equation for u* but using eta_old at old time
    // (we can't really trust the vel we have so far in this step to define eta at new time)
    if (m_diff_type == DiffusionType::Crank_Nicolson)
    {
        diffusion_op->diffuse_velocity(vel   , density, eta_old, 0.5*dt);
        if (advect_tracer)
            diffusion_op->diffuse_scalar  (tracer, density, eta_tracer,    0.5*dt);
    }
    else if (m_diff_type == DiffusionType::Implicit)
    {
        diffusion_op->diffuse_velocity(vel   , density, eta_old, dt);
        if (advect_tracer)
            diffusion_op->diffuse_scalar  (tracer, density, eta_tracer,    dt);
    }

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
    BL_PROFILE("incflo::ApplyCorrector");

    // We use the new time value for things computed on the "*" state
    Real new_time = cur_time + dt;

    if(incflo_verbose > 2)
    {
        amrex::Print() << "Before corrector step:" << std::endl;
        PrintMaxValues(new_time);
    }

    // This fills the eta array (if non-Newtonian, then using strain-rate of velocity at time "new_time",
    //                           which is currently u*)
    // We need this eta whether explicit, implicit or Crank-Nicolson
    ComputeViscosity(eta, new_time);

    // Adds additional viscosity if LES is on
    add_eddy_viscosity(eta, eta_tracer, cur_time);
    
    // fixme move this somewhere else if we want to monitor it
    amrex::Print() << "eta norm " << Norm(eta,0,0,0) << std::endl;
    amrex::Print() << "eta tracer norm " << Norm(eta_tracer,0,0,0) << std::endl;

    
    // Compute explicit diffusion if used -- note that even though we call this "explicit",
    //   the diffusion term does end up being time-centered so formally second-order
    //   Now divtau is the diffusion term computed from u*
    if (m_diff_type == DiffusionType::Explicit)
    {
       int extrap_dir_bcs = 0;
       incflo_set_velocity_bcs (new_time, vel, extrap_dir_bcs);
       diffusion_op->ComputeDivTau(divtau, vel   , density, eta);
       diffusion_op->ComputeLapS  (laps,   tracer, density, mu_s);
    } else {
       for (int lev = 0; lev <= finest_level; lev++)
          divtau[lev]->setVal(0.);
    }

    // Compute the explicit advective terms R_u^* and R_s^*
    incflo_compute_convective_term( conv_u, conv_r, conv_t, vel, density, tracer, new_time );

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

           if (m_diff_type == DiffusionType::Explicit)
           {
               MultiFab::Saxpy(*tracer[lev], dt / 2.0, *laps_old[lev], 0, 0, ntrac, 0);
               MultiFab::Saxpy(*tracer[lev], dt / 2.0,     *laps[lev], 0, 0, ntrac, 0);
           } else if (m_diff_type == DiffusionType::Crank_Nicolson)
               MultiFab::Saxpy(*tracer[lev], dt / 2.0, *laps_old[lev], 0, 0, ntrac, 0);
        }

        // Add the viscous terms         
        if (m_diff_type == DiffusionType::Explicit)
        {
            MultiFab::Saxpy(*vel[lev], dt / 2.0, *divtau_old[lev], 0, 0, AMREX_SPACEDIM, 0);
            MultiFab::Saxpy(*vel[lev], dt / 2.0,     *divtau[lev], 0, 0, AMREX_SPACEDIM, 0);
        } else if (m_diff_type == DiffusionType::Crank_Nicolson)
            MultiFab::Saxpy(*vel[lev], dt / 2.0, *divtau_old[lev], 0, 0, AMREX_SPACEDIM, 0);

        // Add gravitational forces
        if (!use_boussinesq)
        {
           for(int dir = 0; dir < AMREX_SPACEDIM; dir++)
            (*vel[lev]).plus(dt * gravity[dir], dir, 1, 0);
        }

        add_abl_source_terms(*vel[lev],*tracer[lev],vel[lev]->nGrow());
        
        
        // Convert velocities to momenta
        for(int dir = 0; dir < AMREX_SPACEDIM; dir++)
        {
            MultiFab::Multiply(*vel[lev], (*density[lev]), 0, dir, 1, vel[lev]->nGrow());
        }

        // Add (-dt grad p to momenta)
        MultiFab::Saxpy(*vel[lev], -dt, *gp[lev], 0, 0, AMREX_SPACEDIM, vel[lev]->nGrow());

        // Add (-dt grad p0 to momenta if not Boussinesq)
        if (!use_boussinesq)
            for(int dir = 0; dir < AMREX_SPACEDIM; dir++)
                (*vel[lev]).plus(-dt * gp0[dir], dir, 1, 0);

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
    if (m_diff_type == DiffusionType::Crank_Nicolson)
    {
       diffusion_op->diffuse_velocity(vel   , density, eta,  0.5*dt);
       diffusion_op->diffuse_scalar  (tracer, density, eta_tracer, 0.5*dt);
    }
    else if (m_diff_type == DiffusionType::Implicit)
    {
       diffusion_op->diffuse_velocity(vel   , density, eta,  dt);
       diffusion_op->diffuse_scalar  (tracer, density, eta_tracer, dt);
    }

    // Project velocity field, update pressure
    ApplyProjection(new_time, dt);

    // Fill velocity BCs again
    incflo_set_velocity_bcs(new_time, vel, 0);
}
