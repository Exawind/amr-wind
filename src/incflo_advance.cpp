#include <incflo.H>

using namespace amrex;

void incflo::Advance()
{
    BL_PROFILE("incflo::Advance");

    // Start timing current time step
    Real strt_step = ParallelDescriptor::second();

    if (m_verbose > 0)
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
        m_t_old[lev] = m_cur_time;
        m_t_new[lev] = m_cur_time + m_dt;
    }

    if (m_verbose > 0)
    {
        amrex::Print() << "\nStep " << m_nstep + 1
                       << ": from old_time " << m_cur_time
                       << " to new time " << m_cur_time + m_dt
                       << " with dt = " << m_dt << ".\n" << std::endl;
    }

    copy_from_new_to_old_velocity();
    copy_from_new_to_old_density();
    copy_from_new_to_old_tracer();

    int ng = nghost_state();
    for (int lev = 0; lev <= finest_level; ++lev) {
        fillpatch_velocity(lev, m_t_old[lev], m_leveldata[lev]->velocity_o, ng);
        fillpatch_density(lev, m_t_old[lev], m_leveldata[lev]->density_o, ng);
        if (m_advect_tracer) {
            fillpatch_tracer(lev, m_t_old[lev], m_leveldata[lev]->tracer_o, ng);
        }
    }

    ApplyPredictor();

    if (!m_use_godunov) {
        for (int lev = 0; lev <= finest_level; ++lev) {
            fillpatch_velocity(lev, m_t_new[lev], m_leveldata[lev]->velocity, ng);
            fillpatch_density(lev, m_t_new[lev], m_leveldata[lev]->density, ng);
            if (m_advect_tracer) {
                fillpatch_tracer(lev, m_t_new[lev], m_leveldata[lev]->tracer, ng);
            }
        }

        ApplyCorrector();
    }

    if (m_verbose > 1)
    {
        amrex::Print() << "End of time step: " << std::endl;
#if 0
        // xxxxx
        PrintMaxValues(m_cur_time + dt);
        if(m_probtype%10 == 3 or m_probtype == 5)
        {
            ComputeDrag();
            amrex::Print() << "Drag force = " << (*drag[0]).sum(0, false) << std::endl;
        }
#endif
    }

#if 0
    if (m_test_tracer_conservation) {
        amrex::Print() << "Sum tracer volume wgt = " << m_cur_time+dt << "   " << volWgtSum(0,*tracer[0],0) << std::endl;
    }
#endif

    // Stop timing current time step
    Real end_step = ParallelDescriptor::second() - strt_step;
    ParallelDescriptor::ReduceRealMax(end_step, ParallelDescriptor::IOProcessorNumber());
    if (m_verbose > 0)
    {
        amrex::Print() << "Time per step " << end_step << std::endl;
    }
}

//
// Apply predictor:
//
//  1. Use u = vel_old to compute
//
//      conv_u  = - u grad u
//      conv_r  = - div( u rho  )
//      conv_t  = - div( u trac )
//      eta_old     = visosity at m_cur_time
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
    Real new_time = m_cur_time + m_dt;

    if (m_verbose > 2)
    {
        amrex::Print() << "Before predictor step:" << std::endl;
        PrintMaxValues(new_time);
    }

    compute_forces(get_vel_forces(), get_tra_forces(),
                   get_density_old_const(), get_tracer_old_const());

    if (need_divtau()) {
        get_diffusion_tensor_op()->compute_divtau(get_divtau_old(),
                                                  get_velocity_old_const(),
                                                  get_density_old_const(),
                                                  m_cur_time);
        for (int lev = 0; lev <= finest_level; ++lev) {
            MultiFab::Add(m_leveldata[lev]->vel_forces,
                          m_leveldata[lev]->divtau_o, 0, 0, AMREX_SPACEDIM, 0);
        }

        if (m_advect_tracer) {
            get_diffusion_scalar_op()->compute_laps(get_laps_old(),
                                                    get_tracer_old_const(),
                                                    get_density_old_const(),
                                                    m_cur_time);
            for (int lev = 0; lev <= finest_level; ++lev) {
                MultiFab::Add(m_leveldata[lev]->tra_forces,
                              m_leveldata[lev]->laps_o, 0, 0, m_ntrac, 0);
            }
        }
    }

    if (m_use_godunov) {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(Geom(0).isAllPeriodic() and finest_level == 0,
                                         "TODO: fillpatch_forces");
        m_leveldata[0]->vel_forces.FillBoundary(Geom(0).periodicity());
    } else {
        // forces are not used in compute_convective_term
    }

    // if ( m_use_godunov) Compute the explicit advective terms R_u^(n+1/2), R_s^(n+1/2) and R_t^(n+1/2)
    // if (!m_use_godunov) Compute the explicit advective terms R_u^n      , R_s^n       and R_t^n
    compute_convective_term(get_conv_velocity_old(), get_conv_density_old(), get_conv_tracer_old(),
                            get_velocity_old_const(), get_density_old_const(), get_tracer_old_const(),
                            get_vel_forces_const(), get_tra_forces_const(),
                            m_cur_time);

    // Define local variables for lambda to capture.
    Real l_dt = m_dt;
    bool l_constant_density = m_constant_density;
    int l_ntrac = (m_advect_tracer) ? m_ntrac : 0;
    for (int lev = 0; lev <= finest_level; lev++)
    {
        auto& ld = *m_leveldata[lev];
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
            Array4<Real const> const& vel_f = ld.vel_forces.const_array(mfi);
            Array4<Real const> const& tra_f = ld.tra_forces.const_array(mfi);
            // if need_divtau()==true, the forces have already included diffusion terms
            if (need_divtau() and m_diff_type != DiffusionType::Explicit) {
                Array4<Real const> const& divtau = ld.divtau_o.const_array(mfi);
                Array4<Real const> laps;
                if (m_advect_tracer) {
                    laps = ld.laps_o.const_array(mfi);
                }
                // It's either implicit or Crank-Nicolson.
                Real s = (m_diff_type == DiffusionType::Implicit) ? -1.0 : -0.5;
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    vel(i,j,k,0) += l_dt*(dvdt(i,j,k,0)+vel_f(i,j,k,0)+s*divtau(i,j,k,0));
                    vel(i,j,k,1) += l_dt*(dvdt(i,j,k,1)+vel_f(i,j,k,1)+s*divtau(i,j,k,1));
                    vel(i,j,k,2) += l_dt*(dvdt(i,j,k,2)+vel_f(i,j,k,2)+s*divtau(i,j,k,2));

                    if (!l_constant_density) {
                        rho(i,j,k) += l_dt * drdt(i,j,k);
                    }

                    for (int n = 0; n < l_ntrac; ++n) {
                        tra(i,j,k,n) += l_dt * (dtdt(i,j,k,n)+tra_f(i,j,k,n)+s*laps(i,j,k,n));
                    }
                });
            } else {
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    vel(i,j,k,0) += l_dt*(dvdt(i,j,k,0)+vel_f(i,j,k,0));
                    vel(i,j,k,1) += l_dt*(dvdt(i,j,k,1)+vel_f(i,j,k,1));
                    vel(i,j,k,2) += l_dt*(dvdt(i,j,k,2)+vel_f(i,j,k,2));

                    if (!l_constant_density) {
                        rho(i,j,k) += l_dt * drdt(i,j,k);
                    }

                    for (int n = 0; n < l_ntrac; ++n) {
                        tra(i,j,k,n) += l_dt * (dtdt(i,j,k,n)+tra_f(i,j,k,n));
                    }
                });
            }
        }
    }

    // Solve diffusion equation for u* but using eta_old at old time
    if (m_diff_type == DiffusionType::Crank_Nicolson || m_diff_type == DiffusionType::Implicit)
    {
        const int ng_diffusion = 1;
        for (int lev = 0; lev <= finest_level; ++lev) {
            fillphysbc_velocity(lev, new_time, m_leveldata[lev]->velocity, ng_diffusion);
            if (m_advect_tracer) {
                fillphysbc_tracer(lev, new_time, m_leveldata[lev]->tracer, ng_diffusion);
            }
        }

        Real dt_diff = (m_diff_type == DiffusionType::Implicit) ? m_dt : 0.5*m_dt;
        get_diffusion_tensor_op()->diffuse_velocity(get_velocity_new(),
                                                    get_density_new(), m_cur_time, dt_diff);
        if (m_advect_tracer) {
            get_diffusion_scalar_op()->diffuse_scalar(get_tracer_new(),
                                                      get_density_new(), m_cur_time, dt_diff);
        }
    }

    // **********************************************************************************************
    // 
    // Project velocity field, update pressure
    // 
    // **********************************************************************************************
    ApplyProjection(new_time, m_dt, incremental_projection);
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
    Real new_time = m_cur_time + m_dt;

    if (m_verbose > 2)
    {
        amrex::Print() << "Before corrector step:" << std::endl;
        PrintMaxValues(new_time);
    }

    // This fills the eta array (if non-Newtonian, then using strain-rate of velocity at time "new_time",
    //                           which is currently u*)
    // We need this eta whether explicit, implicit or Crank-Nicolson
    if (m_fluid_model != "newtonian") {
        amrex::Abort("non-Newtonian: TODO");
//        ComputeViscosity(eta, new_time);
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
    // We only reach the corrector if !m_use_godunov which means we don't use the forces
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
                            get_velocity_new_const(), get_density_new_const(), get_tracer_new_const(),
                            get_vel_forces_const(), get_tra_forces_const(),
                            new_time);

    compute_forces(get_vel_forces(), get_tra_forces(),
                   get_density_new_const(), get_tracer_new_const());

    // Define local variables for lambda to capture.
    Real l_dt = m_dt;
    bool l_constant_density = m_constant_density;
    int l_ntrac = (m_advect_tracer) ? m_ntrac : 0;
    for (int lev = 0; lev <= finest_level; ++lev)
    {
        auto& ld = *m_leveldata[lev];

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(ld.velocity,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            Box const& bx = mfi.tilebox();
            Array4<Real> const& vel = ld.velocity.array(mfi);
            Array4<Real> const& rho = ld.density.array(mfi);
            Array4<Real> const& tra = ld.tracer.array(mfi);
            Array4<Real const> const& vel_o = ld.velocity_o.const_array(mfi);
            Array4<Real const> const& rho_o = ld.density_o.const_array(mfi);
            Array4<Real const> const& tra_o = ld.tracer_o.const_array(mfi);
            Array4<Real const> const& dvdt = ld.conv_velocity.const_array(mfi);
            Array4<Real const> const& drdt = ld.conv_density.const_array(mfi);
            Array4<Real const> const& dtdt = ld.conv_tracer.const_array(mfi);
            Array4<Real const> const& dvdt_o = ld.conv_velocity_o.const_array(mfi);
            Array4<Real const> const& drdt_o = ld.conv_density_o.const_array(mfi);
            Array4<Real const> const& dtdt_o = ld.conv_tracer_o.const_array(mfi);
            Array4<Real const> const& vel_f = ld.vel_forces.const_array(mfi);
            Array4<Real const> const& tra_f = ld.tra_forces.const_array(mfi);

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                    vel(i,j,k,idim) = vel_o(i,j,k,idim) + l_dt*
                        (0.5*(dvdt_o(i,j,k,idim)+dvdt(i,j,k,idim))+vel_f(i,j,k,idim));
                }

                if (!l_constant_density) {
                    rho(i,j,k) = rho_o(i,j,k) + l_dt * 0.5*(drdt(i,j,k)+drdt_o(i,j,k));
                }

                for (int n = 0; n < l_ntrac; ++n) {
                    tra(i,j,k,n) = tra_o(i,j,k,n) + l_dt *
                        (0.5*(dtdt(i,j,k,n)+dtdt_o(i,j,k,n)+tra_f(i,j,k,n)));
                }
            });
        }
    }

    // **********************************************************************************************
    // 
    // Solve diffusion equation for u* at t^{n+1} but using eta at predicted new time
    // 
    // **********************************************************************************************

    if (m_diff_type == DiffusionType::Crank_Nicolson)
    {
        amrex::Abort("TODO: Explicit or Crank_Nicolson diffustion");
        // diffusion_op->diffuse_velocity(vel   , density, eta,  0.5*dt);
        // diffusion_op->diffuse_scalar  (tracer, density, mu_s, 0.5*dt);
    }
    else if (m_diff_type == DiffusionType::Implicit)
    {
        get_diffusion_tensor_op()->diffuse_velocity(get_velocity_new(),
                                                    get_density_new(), new_time, m_dt);
        if (m_advect_tracer) {
            get_diffusion_scalar_op()->diffuse_scalar(get_tracer_new(),
                                                      get_density_new(), new_time, m_dt);
        }
    }

    // **********************************************************************************************
    // 
    // Project velocity field, update pressure
    bool incremental = false;
    ApplyProjection(new_time, m_dt, incremental);
}

void incflo::compute_forces (Vector<MultiFab*> const& vel_forces,
                             Vector<MultiFab*> const& tra_forces,
                             Vector<MultiFab const*> const& density,
                             Vector<MultiFab const*> const& tracer)
{
    GpuArray<Real,3> l_gravity{m_gravity[0],m_gravity[1],m_gravity[2]};
    GpuArray<Real,3> l_gp0{m_gp0[0], m_gp0[1], m_gp0[2]};

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (int lev = 0; lev <= finest_level; ++lev) {
        for (MFIter mfi(*vel_forces[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            Box const& bx = mfi.tilebox();
            Array4<Real> const& vel_f = vel_forces[lev]->array(mfi);
            Array4<Real const> const& rho = density[lev]->const_array(mfi);
            Array4<Real const> const& gradp = m_leveldata[lev]->gp.const_array(mfi);

            if (m_use_boussinesq) {
                // This uses a Boussinesq approximation where the buoyancy depends on
                //      first tracer rather than density
                Array4<Real const> const& tra = tracer[lev]->const_array(mfi);
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    Real rhoinv = 1.0/rho(i,j,k);
                    Real ft = tra(i,j,k);
                    vel_f(i,j,k,0) = -gradp(i,j,k,0)*rhoinv + l_gravity[0] * ft;
                    vel_f(i,j,k,1) = -gradp(i,j,k,1)*rhoinv + l_gravity[1] * ft;
                    vel_f(i,j,k,2) = -gradp(i,j,k,2)*rhoinv + l_gravity[2] * ft;
                });
            } else {
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    Real rhoinv = 1.0/rho(i,j,k);
                    vel_f(i,j,k,0) = -(gradp(i,j,k,0)+l_gp0[0])*rhoinv + l_gravity[0];
                    vel_f(i,j,k,1) = -(gradp(i,j,k,1)+l_gp0[1])*rhoinv + l_gravity[1];
                    vel_f(i,j,k,2) = -(gradp(i,j,k,2)+l_gp0[2])*rhoinv + l_gravity[2];
                });
            }
        }
    }

    // For now we don't have any external forces on the scalars
    for (int lev = 0; lev <= finest_level; ++lev) {
        tra_forces[lev]->setVal(0.0);
    }
}
