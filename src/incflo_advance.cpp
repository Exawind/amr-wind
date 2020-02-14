#include <incflo.H>
#include "PlaneAveraging.H"
#include <cmath>

using namespace amrex;

void incflo::Advance()
{
    BL_PROFILE("incflo::Advance");

    // Start timing current time step
    Real strt_step = ParallelDescriptor::second();

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

   if(m_plane_averaging){
       const int axis=2;
       PlaneAveraging pa(Geom(), get_velocity_new(), get_tracer_new(), axis);

       Real vx = pa.line_velocity_xdir(geom[0].ProbLoArray()[axis]);
       Real vy = pa.line_velocity_ydir(geom[0].ProbLoArray()[axis]);

       m_velocity_mean_ground = std::sqrt(vx*vx + vy*vy);

       vx = pa.line_velocity_xdir(m_log_law_sampling_height);
       vy = pa.line_velocity_ydir(m_log_law_sampling_height);

       m_velocity_mean_loglaw = std::sqrt(vx*vx + vy*vy);

       m_vx_mean_forcing = pa.line_velocity_xdir(m_abl_forcing_height);
       m_vy_mean_forcing = pa.line_velocity_ydir(m_abl_forcing_height);

       if(m_line_plot_int > 0 and m_nstep % m_line_plot_int == 0)
       {
           pa.plot_line_text("line_plot.txt", m_nstep, m_cur_time);
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

    if (m_verbose > 2)
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

    Vector<MultiFab> vel_forces, tra_forces;
    Vector<MultiFab> vel_eta, tra_eta;
    for (int lev = 0; lev <= finest_level; ++lev) {
        vel_forces.emplace_back(grids[lev], dmap[lev], AMREX_SPACEDIM, nghost_force(),
                                MFInfo(), Factory(lev));

        if (m_advect_tracer) {
            tra_forces.emplace_back(grids[lev], dmap[lev], m_ntrac, nghost_force(),
                                    MFInfo(), Factory(lev));
        }
        vel_eta.emplace_back(grids[lev], dmap[lev], 1, 1, MFInfo(), Factory(lev));
        if (m_advect_tracer) {
            tra_eta.emplace_back(grids[lev], dmap[lev], m_ntrac, 1, MFInfo(), Factory(lev));
        }
    }

    compute_forces(GetVecOfPtrs(vel_forces), GetVecOfPtrs(tra_forces),
                   get_velocity_old_const(), get_density_old_const(), get_tracer_old_const(),m_dt);

    compute_viscosity(GetVecOfPtrs(vel_eta), GetVecOfPtrs(tra_eta),
                      get_density_old_const(), get_velocity_old_const(), get_tracer_old_const(),
                      m_cur_time, 1);

    if (need_divtau()) {
        get_diffusion_tensor_op()->compute_divtau(get_divtau_old(),
                                                  get_velocity_old_const(),
                                                  get_density_old_const(),
                                                  GetVecOfConstPtrs(vel_eta),
                                                  m_cur_time);
        for (int lev = 0; lev <= finest_level; ++lev) {
            MultiFab::Add(vel_forces[lev], m_leveldata[lev]->divtau_o, 0, 0, AMREX_SPACEDIM, 0);
        }

        if (m_advect_tracer) {
            get_diffusion_scalar_op()->compute_laps(get_laps_old(),
                                                    get_tracer_old_const(),
                                                    get_density_old_const(),
                                                    GetVecOfConstPtrs(tra_eta),
                                                    m_cur_time);
            for (int lev = 0; lev <= finest_level; ++lev) {
                MultiFab::Add(tra_forces[lev], m_leveldata[lev]->laps_o, 0, 0, m_ntrac, 0);
            }
        }
    }

    if (m_use_godunov and nghost_force() > 0) {
        fillpatch_force(m_cur_time, GetVecOfPtrs(vel_forces), nghost_force());
        if (m_advect_tracer) {
            fillpatch_force(m_cur_time, GetVecOfPtrs(tra_forces), nghost_force());
        }
    }

    // if ( m_use_godunov) Compute the explicit advective terms R_u^(n+1/2), R_s^(n+1/2) and R_t^(n+1/2)
    // if (!m_use_godunov) Compute the explicit advective terms R_u^n      , R_s^n       and R_t^n
    compute_convective_term(get_conv_velocity_old(), get_conv_density_old(), get_conv_tracer_old(),
                            get_velocity_old_const(), get_density_old_const(), get_tracer_old_const(),
                            GetVecOfConstPtrs(vel_forces), GetVecOfConstPtrs(tra_forces),
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
            Array4<Real const> const& vel_f = vel_forces[lev].const_array(mfi);
            Array4<Real const> const& tra_f = (l_ntrac > 0) ? tra_forces[lev].const_array(mfi)
                                                            : Array4<Real const>{};
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
                                                    get_density_new(),
                                                    GetVecOfConstPtrs(vel_eta),
                                                    m_cur_time, dt_diff);
        if (m_advect_tracer) {
            get_diffusion_scalar_op()->diffuse_scalar(get_tracer_new(),
                                                      get_density_new(),
                                                      GetVecOfConstPtrs(tra_eta),
                                                      m_cur_time, dt_diff);
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

    // **********************************************************************************************
    // 
    // We only reach the corrector if !m_use_godunov which means we don't use the forces
    //    in constructing the advection term
    // 
    // **********************************************************************************************
    Vector<MultiFab> vel_forces, tra_forces;
    Vector<MultiFab> vel_eta, tra_eta;
    for (int lev = 0; lev <= finest_level; ++lev) {
        vel_forces.emplace_back(grids[lev], dmap[lev], AMREX_SPACEDIM, nghost_force(),
                                MFInfo(), Factory(lev));
        if (m_advect_tracer) {
            tra_forces.emplace_back(grids[lev], dmap[lev], m_ntrac, nghost_force(),
                                    MFInfo(), Factory(lev));
        }
        vel_eta.emplace_back(grids[lev], dmap[lev], 1, 1, MFInfo(), Factory(lev));
        if (m_advect_tracer) {
            tra_eta.emplace_back(grids[lev], dmap[lev], m_ntrac, 1, MFInfo(), Factory(lev));
        }
    }

    // **********************************************************************************************
    // 
    // Compute convective / conservative update
    // 
    // **********************************************************************************************

    compute_convective_term(get_conv_velocity_new(), get_conv_density_new(), get_conv_tracer_new(),
                            get_velocity_new_const(), get_density_new_const(), get_tracer_new_const(),
                            {}, {}, new_time);

    compute_forces(GetVecOfPtrs(vel_forces), GetVecOfPtrs(tra_forces),
                   get_velocity_new_const(), get_density_new_const(), get_tracer_new_const(),m_dt);

    compute_viscosity(GetVecOfPtrs(vel_eta), GetVecOfPtrs(tra_eta),
                      get_density_new_const(), get_velocity_new_const(), get_tracer_new_const(),
                      new_time, 1);

    if (m_diff_type == DiffusionType::Explicit) {
        get_diffusion_tensor_op()->compute_divtau(get_divtau_new(),
                                                  get_velocity_new_const(),
                                                  get_density_new_const(),
                                                  GetVecOfConstPtrs(vel_eta),
                                                  m_cur_time);
        if (m_advect_tracer) {
            get_diffusion_scalar_op()->compute_laps(get_laps_new(),
                                                    get_tracer_new_const(),
                                                    get_density_new_const(),
                                                    GetVecOfConstPtrs(tra_eta),
                                                    m_cur_time);
        }
    }

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
            Array4<Real const> const& vel_f = vel_forces[lev].const_array(mfi);
            Array4<Real const> const& tra_f = (l_ntrac > 0) ? tra_forces[lev].const_array(mfi)
                                                            : Array4<Real const>{};

            if (m_diff_type == DiffusionType::Explicit) {
                Array4<Real const> const& divtau_o = ld.divtau_o.const_array(mfi);
                Array4<Real const> const& divtau   = ld.divtau.const_array(mfi);
                Array4<Real const> const& laps_o = (l_ntrac > 0) ? ld.laps_o.const_array(mfi)
                                                                 : Array4<Real const>{};
                Array4<Real const> const& laps   = (l_ntrac > 0) ? ld.laps.const_array(mfi)
                                                                 : Array4<Real const>{};
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                        vel(i,j,k,idim) = vel_o(i,j,k,idim) + l_dt*
                            (0.5*(dvdt_o(i,j,k,idim)+dvdt(i,j,k,idim)
                                +divtau_o(i,j,k,idim)+divtau(i,j,k,idim))
                             +vel_f(i,j,k,idim));
                    }

                    if (!l_constant_density) {
                        rho(i,j,k) = rho_o(i,j,k) + l_dt * 0.5*(drdt(i,j,k)+drdt_o(i,j,k));
                    }

                    for (int n = 0; n < l_ntrac; ++n) {
                        tra(i,j,k,n) = tra_o(i,j,k,n) + l_dt *
                            (0.5*(dtdt(i,j,k,n)+dtdt_o(i,j,k,n)+laps_o(i,j,k,n)+laps(i,j,k,n))
                             +tra_f(i,j,k,n));
                    }
                });
            } else if (m_diff_type == DiffusionType::Crank_Nicolson) {
                Array4<Real const> const& divtau_o = ld.divtau_o.const_array(mfi);
                Array4<Real const> const& laps_o = (l_ntrac > 0) ? ld.laps_o.const_array(mfi)
                                                                 : Array4<Real const>{};
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                        vel(i,j,k,idim) = vel_o(i,j,k,idim) + l_dt*
                            (0.5*(dvdt_o(i,j,k,idim)+dvdt(i,j,k,idim))+vel_f(i,j,k,idim) 
                                +divtau_o(i,j,k,idim));
                    }

                    if (!l_constant_density) {
                        rho(i,j,k) = rho_o(i,j,k) + l_dt * 0.5*(drdt(i,j,k)+drdt_o(i,j,k));
                    }

                    for (int n = 0; n < l_ntrac; ++n) {
                        tra(i,j,k,n) = tra_o(i,j,k,n) + l_dt *
                            (0.5*(dtdt(i,j,k,n)+dtdt_o(i,j,k,n))+tra_f(i,j,k,n)
                                +laps_o(i,j,k,n));
                    }
                });
            } else {
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
                            (0.5*(dtdt(i,j,k,n)+dtdt_o(i,j,k,n))+tra_f(i,j,k,n));
                    }
                });
            }
        }
    }

    // **********************************************************************************************
    // 
    // Solve diffusion equation for u* at t^{n+1} but using eta at predicted new time
    // 
    // **********************************************************************************************

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
                                                    get_density_new(),
                                                    GetVecOfConstPtrs(vel_eta),
                                                    new_time, dt_diff);
        if (m_advect_tracer) {
            get_diffusion_scalar_op()->diffuse_scalar(get_tracer_new(),
                                                      get_density_new(),
                                                      GetVecOfConstPtrs(tra_eta),
                                                      new_time, dt_diff);
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
                             Vector<MultiFab const*> const& velocity,
                             Vector<MultiFab const*> const& density,
                             Vector<MultiFab const*> const& tracer,
                             Real dt)
{
    // For now we don't have any external forces on the scalars
    if (m_advect_tracer) {
        for (int lev = 0; lev <= finest_level; ++lev) {
            tra_forces[lev]->setVal(0.0);
        }
    }

    GpuArray<Real,3> l_gravity{m_gravity[0],m_gravity[1],m_gravity[2]};
    GpuArray<Real,3> l_gp0{m_gp0[0], m_gp0[1], m_gp0[2]};
    
    GpuArray<Real,3> east = {m_east[0],m_east[1],m_east[2]};
    GpuArray<Real,3> north = {m_north[0],m_north[1],m_north[2]};
    GpuArray<Real,3> up = {m_up[0],m_up[1],m_up[2]};
     
    const Real sinphi = m_sinphi;
    const Real cosphi = m_cosphi;
    const Real corfac = m_corfac;

    const Real u = m_ic_u;
    const Real v = m_ic_v;
    
    const Real umean = m_vx_mean_forcing;//fixme get rid of this global storage
    const Real vmean = m_vy_mean_forcing;//fixme get rid of this global storage
    
    const Real T0 = m_temperature_values[0];
    const Real thermalExpansionCoeff = m_thermalExpansionCoeff;
    AMREX_ALWAYS_ASSERT(thermalExpansionCoeff>0.0);
    
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
                    Real ft = thermalExpansionCoeff*(T0-tra(i,j,k));
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
            
            if(m_use_coriolis){
                
                Array4<Real const> const& vel = velocity[lev]->const_array(mfi);
                
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    
                    const Real ue = east[0]*vel(i,j,k,0) + east[1]*vel(i,j,k,1) + east[2]*vel(i,j,k,2);
                    const Real un = north[0]*vel(i,j,k,0) + north[1]*vel(i,j,k,1) + north[2]*vel(i,j,k,2);
                    const Real uu = up[0]*vel(i,j,k,0) + up[1]*vel(i,j,k,1) + up[2]*vel(i,j,k,2);

                    const Real ae = +corfac*(un*sinphi - uu*cosphi);
                    const Real an = -corfac*ue*sinphi;
                    const Real au = +corfac*ue*cosphi;

                    const Real ax = ae*east[0] + an*north[0] + au*up[0];
                    const Real ay = ae*east[1] + an*north[1] + au*up[1];
                    const Real az = ae*east[2] + an*north[2] + au*up[2];//fixme do we turn the z component off?
                    
                    vel_f(i,j,k,0) += ax;
                    vel_f(i,j,k,1) += ay;
                    vel_f(i,j,k,2) += az;
                    
                });
            }
            
            if(m_use_abl_forcing){
                const Real dudt = (u-umean)/dt;// fixme make sure this is the correct dt and not dt/2
                const Real dvdt = (v-vmean)/dt;
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    vel_f(i,j,k,0) += dudt;
                    vel_f(i,j,k,1) += dvdt;
                    // vel_f(i,j,k,2) -= 0.0; // assuming periodic in x,y-dir so do not drive flow in z-dir
                });
            }
            
        }
    }
}
