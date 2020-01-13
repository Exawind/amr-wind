#include <incflo.H>

using namespace amrex;

void incflo::Advance()
{
    BL_PROFILE("incflo::Advance");

    // Start timing current time step
    Real strt_step = ParallelDescriptor::second();

    if(incflo_verbose > 0)
    {
        amrex::Print() << "\n ============   NEW TIME STEP   ============ \n";
    }

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

    copy_from_new_to_old_velocity();
    copy_from_new_to_old_density();
    copy_from_new_to_old_tracer();

    int ng = 2;  // This might change for Godunov.
#ifdef AMREX_USE_EB
    if (!EBFactory(0).isAllRegular()) ng = 4;
#endif
    for (int lev = 0; lev <= finest_level; ++lev) {
        fillpatch_velocity(lev, t_old[lev], m_leveldata[lev]->velocity_o, ng);
        fillpatch_density(lev, t_old[lev], m_leveldata[lev]->density_o, ng);
        if (advect_tracer) {
            fillpatch_tracer(lev, t_old[lev], m_leveldata[lev]->tracer_o, ng);
        }
    }

    ApplyPredictor();

    for (int lev = 0; lev <= finest_level; ++lev) {
        fillpatch_velocity(lev, t_new[lev], m_leveldata[lev]->velocity, ng);
        fillpatch_density(lev, t_new[lev], m_leveldata[lev]->density, ng);
        if (advect_tracer) {
            fillpatch_tracer(lev, t_new[lev], m_leveldata[lev]->tracer, ng);
        }
    }

    if (!use_godunov) {
        ApplyCorrector();
    }

    amrex::EB_set_covered(m_leveldata[0]->velocity, 0.0);
    amrex::EB_set_covered(m_leveldata[0]->density, 0.0);
    amrex::EB_set_covered(m_leveldata[0]->tracer, 0.0);
    amrex::EB_set_covered(m_leveldata[0]->p, 0.0);
    amrex::VisMF::Write(m_leveldata[0]->velocity, "vel");
    amrex::VisMF::Write(m_leveldata[0]->density, "rho");
    amrex::VisMF::Write(m_leveldata[0]->tracer, "tra");
    amrex::VisMF::Write(m_leveldata[0]->p, "p");

    amrex::Abort("xxxxx in Advance(): after ApplyCorrector");

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

    if (test_tracer_conservation)
       amrex::Print() << "Sum tracer volume wgt = " << cur_time+dt << "   " << volWgtSum(0,*tracer[0],0) << std::endl;

    // Stop timing current time step
    Real end_step = ParallelDescriptor::second() - strt_step;
    ParallelDescriptor::ReduceRealMax(end_step, ParallelDescriptor::IOProcessorNumber());
    if(incflo_verbose > 0)
    {
        amrex::Print() << "Time per step " << end_step << std::endl;
    }

    amrex::Abort("xxxxx in Advance(): end");
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
// It is assumed that the ghost cels of the old data have been filled and
// the old and new data are the same in valid region.
//
void incflo::ApplyPredictor (bool incremental_projection)
{
    BL_PROFILE("incflo::ApplyPredictor");

    // We use the new time value for things computed on the "*" state
    Real new_time = cur_time + dt;

    if(incflo_verbose > 2)
    {
        amrex::Print() << "Before predictor step:" << std::endl;
        PrintMaxValues(new_time);
    }

    // if ( use_godunov) Compute the explicit advective terms R_u^(n+1/2), R_s^(n+1/2) and R_t^(n+1/2)
    // if (!use_godunov) Compute the explicit advective terms R_u^n      , R_s^n       and R_t^n
    compute_convective_term(get_conv_velocity_old(), get_conv_density_old(), get_conv_tracer_old(),
                            get_velocity_old(), get_density_old(), get_tracer_old(), cur_time);

    // This fills the eta_old array (if non-Newtonian, then using strain-rate of velocity at time "cur_time")
    if (fluid_model != "newtonian") {
        amrex::Abort("non-Newtonian: TODO");
        ComputeViscosity(eta_old, cur_time);
    }

    // Compute explicit diffusion if used
    if (m_diff_type == DiffusionType::Explicit ||
        m_diff_type == DiffusionType::Crank_Nicolson ||
        use_godunov)
    {
        amrex::Abort("TODO: Explicit or Crank_Nicolson diffustion");
        // incflo_set_velocity_bcs (cur_time, vel_o);
        // diffusion_op->ComputeDivTau(divtau_old,    vel_o, density_o, eta_old);
        // diffusion_op->ComputeLapS  (  laps_old, tracer_o, density_o, mu_s);
    } else {
#if 0
        // xxxxx TODO
       for (int lev = 0; lev <= finest_level; lev++)
       {
          divtau_old[lev]->setVal(0.);
            laps_old[lev]->setVal(0.);
       }
#endif
    }

    // Define local variables for lambda to capture.
    Real l_dt = dt;
    bool l_constant_density = constant_density;
    bool l_use_boussinesq = use_boussinesq;
    int l_ntrac = (advect_tracer) ? incflo::ntrac : 0;
    GpuArray<Real,3> l_gravity{gravity[0],gravity[1],gravity[2]};
    GpuArray<Real,3> l_gp0{gp0[0], gp0[1], gp0[2]};

    for (int lev = 0; lev <= finest_level; lev++)
    {
        auto& ld = *m_leveldata[lev];
#if 0
        // xxxxx TODO: conflict
        // Fill forcing term with -gp to start
        MultiFab::Copy (*vel_forces[lev], *gp[lev], 0, 0, 3, 2);
        vel_forces[lev] -> mult(-1.0);

        // For now we don't have any external forces acting on the scalars
        scal_forces[lev] -> setVal(0.);

        // Add (-gp0 to forcing term if not Boussinesq)
        if (!use_boussinesq)
            for(int dir = 0; dir < AMREX_SPACEDIM; dir++)
                (*vel_forces[lev]).plus(-gp0[dir], dir, 1, vel_forces[lev]->nGrow());
#endif

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(ld.velocity,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            Box const& bx = mfi.tilebox();
            Array4<Real> const& vel = ld.velocity.array(mfi);
            Array4<Real> const& rho = ld.density.array(mfi);
            Array4<Real> const& tra = ld.tracer.array(mfi);
            Array4<Real const> const& dvdt = ld.conv_velocity_o.const_array(mfi);
            Array4<Real const> const& drdt = ld.conv_density_o.const_array(mfi);
            Array4<Real const> const& dtdt = ld.conv_tracer_o.const_array(mfi);
            Array4<Real const> const& gradp = ld.gp.const_array(mfi);

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                // Now add the convective term to the density and tracer
                if (!l_constant_density) {
                    rho(i,j,k) += l_dt * drdt(i,j,k);
                }

                if (l_ntrac > 0) {
                    // xxxxx TODO: this is currently inconsistent.
                    //             dtdt is actually for rho*tracer.
                    for (int n = 0; n < l_ntrac; ++n) {
                        tra(i,j,k,n) += l_dt * dtdt(i,j,k,n);
                    }
                }

                // xxxxx TODO add viscous terms for explicit and Crank_Nicolson diffusion type

                Real rhoinv = 1.0/(rho(i,j,k)+1.e-80);
                for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
                    // First add the convective term to the velocity
                    Real vt = dvdt(i,j,k,dir);

                    // Add gravitational forces
                    if (l_use_boussinesq) {
                        // xxxxx TODO: should this use the old or new tracer?
                        vt += tra(i,j,k,0) * l_gravity[dir];
                    } else {
                        vt += l_gravity[dir];
                    }

                    // Add -gradp / rho
                    // xxxxx TODO: should this use the old or new density?
                    vt += -gradp(i,j,k,dir) * rhoinv;

                    if (!l_use_boussinesq) {
                        vt += -l_gp0[dir];
                    }

                    vel(i,j,k,dir) += l_dt * vt;
                }
            });
        }
    }

    if (m_diff_type == DiffusionType::Crank_Nicolson or
        m_diff_type == DiffusionType::Implicit)
    {
        const int ng_diffusion = 1;
        for (int lev = 0; lev <= finest_level; ++lev) {
            fillphysbc_velocity(lev, new_time, m_leveldata[lev]->velocity, ng_diffusion);
            if (advect_tracer) {
                fillphysbc_tracer(lev, new_time, m_leveldata[lev]->tracer, ng_diffusion);
            }
        }
    }

    // Solve diffusion equation for u* but using eta_old at old time
    // (we can't really trust the vel we have so far in this step to define eta at new time)
    // 
    // **********************************************************************************************
    if (m_diff_type == DiffusionType::Crank_Nicolson)
    {
        amrex::Abort("TODO: Crank_Nicolson");
#if 0
        diffusion_op->diffuse_velocity(vel   , density, eta_old, 0.5*dt);
        if (advect_tracer)
            diffusion_op->diffuse_scalar  (tracer, density, mu_s,    0.5*dt);
#endif
    }
    else if (m_diff_type == DiffusionType::Implicit)
    {
        get_diffusion_tensor_op()->diffuse_velocity(get_velocity_new(),
                                                    get_density_new(), cur_time, dt);
        if (advect_tracer) {
            get_diffusion_scalar_op()->diffuse_scalar(get_tracer_new(),
                                                      get_density_new(), cur_time, dt);
        }
    }

    // **********************************************************************************************
    // 
    // Project velocity field, update pressure
    // 
    // **********************************************************************************************
    ApplyProjection(new_time, dt, incremental_projection);
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
    if (fluid_model != "newtonian") {
        amrex::Abort("non-Newtonian: TODO");
        ComputeViscosity(eta, new_time);
    }

    // Compute explicit diffusion if used -- note that even though we call this "explicit",
    //   the diffusion term does end up being time-centered so formally second-order
    //   Now divtau is the diffusion term computed from u*
    if (m_diff_type == DiffusionType::Explicit)
    {
        amrex::Abort("TODO: Explicit");
#if 0
       incflo_set_velocity_bcs (new_time, vel);
       diffusion_op->ComputeDivTau(divtau, vel   , density, eta);
       diffusion_op->ComputeLapS  (laps,   tracer, density, mu_s);
#endif
    } else {
#if 0
       for (int lev = 0; lev <= finest_level; lev++)
          divtau[lev]->setVal(0.);
#endif
    }


    // **********************************************************************************************
    // 
    // We only reach the corrector if !use_godunov which means we don't use the forces
    //    in constructing the advection term
    // 
    // **********************************************************************************************
#if 0
    for (int lev = 0; lev <= finest_level; lev++)
    {
         vel_forces[lev] -> setVal(0.);
        scal_forces[lev] -> setVal(0.);
    }
#endif

    // **********************************************************************************************
    // 
    // Compute convective / conservative update
    // 
    // **********************************************************************************************

    compute_convective_term(get_conv_velocity_new(), get_conv_density_new(), get_conv_tracer_new(),
                            get_velocity_new(), get_density_new(), get_tracer_new(), new_time);


    amrex::EB_set_covered(m_leveldata[0]->conv_velocity, 0.0);
    amrex::EB_set_covered(m_leveldata[0]->conv_density, 0.0);
    amrex::EB_set_covered(m_leveldata[0]->conv_tracer, 0.0);
    VisMF::Write(m_leveldata[0]->conv_velocity, "conv_vel");
    VisMF::Write(m_leveldata[0]->conv_density, "conv_rho");
    VisMF::Write(m_leveldata[0]->conv_tracer, "conv_tra");
    amrex::Abort("xxxxx");

    // **********************************************************************************************
    // 
    // Compute forcing terms to be added explicitly in the update
    // We need the forcing terms to have two ghost cells because of the way they are used 
    //    in the Godunov prediction
    // 
    // **********************************************************************************************
    for (int lev = 0; lev <= finest_level; lev++)
    {
        // Fill forcing term with -gp to start
        MultiFab::Copy (*vel_forces[lev], *gp[lev], 0, 0, 3, 2);
        vel_forces[lev] -> mult(-1.0);

        // Add (-gp0 to forcing term if not Boussinesq)
        if (!use_boussinesq)
            for(int dir = 0; dir < AMREX_SPACEDIM; dir++)
                (*vel_forces[lev]).plus(-gp0[dir], dir, 1, vel_forces[lev]->nGrow());

        // Now divide -(gp+gp0) by density at half-time
        // TODO: CREATE HALF_TIME DENSITY ARRAY
        for(int dir = 0; dir < AMREX_SPACEDIM; dir++)
            MultiFab::Divide(*vel_forces[lev], (*density[lev]), 0, dir, 1, vel_forces[lev]->nGrow());

        // Add gravitational forces
        if (use_boussinesq)
        {
           // TODO: CREATE HALF_TIME TRACER ARRAY
           // This uses a Boussinesq approximation where the buoyancy depends on first tracer
           //      rather than density
           for(int dir = 0; dir < AMREX_SPACEDIM; dir++)
              MultiFab::Saxpy(*vel_forces[lev], gravity[dir], *tracer[lev], 0, dir, 1, 0);
        } else {
           for(int dir = 0; dir < AMREX_SPACEDIM; dir++)
               (*vel_forces[lev]).plus(gravity[dir], dir, 1, 0);
        }

       // Add the diffusive/viscous term to the forcing term 
       if (m_diff_type == DiffusionType::Explicit)
       {
          // Note that even though we call this "explicit",
          //   the diffusion term does end up being time-centered so formally second-order
          MultiFab::Saxpy( *vel_forces[lev],0.5, *divtau_old[lev], 0, 0, AMREX_SPACEDIM, 0);
          MultiFab::Saxpy( *vel_forces[lev],0.5,     *divtau[lev], 0, 0, AMREX_SPACEDIM, 0);

          MultiFab::Saxpy(*scal_forces[lev],0.5,   *laps_old[lev], 0, 0, ntrac         , 0);
          MultiFab::Saxpy(*scal_forces[lev],0.5,       *laps[lev], 0, 0, ntrac         , 0);
       }
       else if (m_diff_type == DiffusionType::Crank_Nicolson)
       {
          MultiFab::Saxpy( *vel_forces[lev],1.0, *divtau_old[lev], 0, 0, AMREX_SPACEDIM, 0);
          MultiFab::Saxpy(*scal_forces[lev],1.0,   *laps_old[lev], 0, 0, ntrac         , 0);
       }

        vel_forces[lev]->FillBoundary(geom[lev].periodicity());
       scal_forces[lev]->FillBoundary(geom[lev].periodicity());
    }

    // **********************************************************************************************
    // 
    // Add advection and other forcing terms to the density, tracer and velocity
    // 
    // **********************************************************************************************
    for(int lev = 0; lev <= finest_level; lev++)
    {
        // Aadd the time-centered advection term to velocity
        MultiFab::LinComb(*vel[lev], 1.0   , *vel_o[lev]     , 0, 0.5*dt, *conv_u[lev], 0, 0, AMREX_SPACEDIM, 0);
        MultiFab::Saxpy  (*vel[lev], 0.5*dt, *conv_u_old[lev], 0, 0, AMREX_SPACEDIM, 0);

        // Add the time-centered forcing term to velocity
        MultiFab::Saxpy(*vel[lev], dt, *vel_forces[lev], 0, 0, AMREX_SPACEDIM, vel_forces[lev]->nGrow());

        if (!constant_density)
        {
            // Add the time-centered advection terms to tracer
            MultiFab::LinComb(*density[lev], 1.0   , *density_o[lev] , 0, 0.5*dt, *conv_r[lev], 0, 0, 1, 0);
            MultiFab::Saxpy  (*density[lev], 0.5*dt, *conv_r_old[lev], 0, 0, 1, 0);
        }

        if (advect_tracer)
        {
            // Add the time-centered advection term to tracer
            MultiFab::LinComb(*tracer[lev], 1.0   , *tracer_o[lev]  , 0, 0.5*dt, *conv_t[lev]    , 0, 0, tracer[lev]->nComp(), 0);
            MultiFab::Saxpy(  *tracer[lev], 0.5*dt, *conv_t_old[lev], 0, 0, tracer[lev]->nComp(), 0);

            // Add the time-centered forcing term to tracer
            MultiFab::Saxpy(*tracer[lev], dt, *scal_forces[lev], 0, 0, ntrac, scal_forces[lev]->nGrow());
        }
    }

    if (!constant_density)
       incflo_set_density_bcs(new_time, density);
    if (advect_tracer)
       incflo_set_tracer_bcs(new_time, tracer);
    incflo_set_velocity_bcs(new_time, vel);

    // **********************************************************************************************
    // 
    // Solve diffusion equation for u* at t^{n+1} but using eta at predicted new time
    // 
    // **********************************************************************************************

    if (m_diff_type == DiffusionType::Crank_Nicolson)
    {
       diffusion_op->diffuse_velocity(vel   , density, eta,  0.5*dt);
       diffusion_op->diffuse_scalar  (tracer, density, mu_s, 0.5*dt);
    }
    else if (m_diff_type == DiffusionType::Implicit)
    {
       diffusion_op->diffuse_velocity(vel   , density, eta,  dt);
       diffusion_op->diffuse_scalar  (tracer, density, mu_s, dt);
    }

    // **********************************************************************************************
    // 
    // Project velocity field, update pressure
    bool incremental = false;
    ApplyProjection(new_time, dt, incremental);

    // Fill velocity BCs again
    incflo_set_velocity_bcs(new_time, vel);
}
