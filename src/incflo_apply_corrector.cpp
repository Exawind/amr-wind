#include <incflo.H>

using namespace amrex;
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

    // *************************************************************************************
    // Allocate space for the MAC velocities
    // *************************************************************************************
    Vector<MultiFab> u_mac(finest_level+1), v_mac(finest_level+1), w_mac(finest_level+1);
    int ngmac = nghost_mac();

    for (int lev = 0; lev <= finest_level; ++lev) {
        u_mac[lev].define(amrex::convert(grids[lev],IntVect::TheDimensionVector(0)), dmap[lev],
                          1, ngmac, MFInfo(), Factory(lev));
        v_mac[lev].define(amrex::convert(grids[lev],IntVect::TheDimensionVector(1)), dmap[lev],
                          1, ngmac, MFInfo(), Factory(lev));
        w_mac[lev].define(amrex::convert(grids[lev],IntVect::TheDimensionVector(2)), dmap[lev],
                          1, ngmac, MFInfo(), Factory(lev));
        if (ngmac > 0) {
            u_mac[lev].setBndry(0.0);
            v_mac[lev].setBndry(0.0);
            w_mac[lev].setBndry(0.0);
        }
    }

    // *************************************************************************************
    // Allocate space for half-time density
    // *************************************************************************************
    Vector<MultiFab> density_nph;
    for (int lev = 0; lev <= finest_level; ++lev) {
        density_nph.emplace_back(grids[lev], dmap[lev], 1, 0, MFInfo(), Factory(lev));
    }

    // **********************************************************************************************
    // We only reach the corrector if !m_use_godunov which means we don't use the forces
    //    in constructing the advection term
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
    // Compute the explicit "new" advective terms R_u^(n+1,*), R_r^(n+1,*) and R_t^(n+1,*)
    // Note that "get_conv_tracer_new" returns div(rho u tracer)
    // *************************************************************************************
    compute_convective_term(get_conv_velocity_new(), get_conv_density_new(), get_conv_tracer_new(),
                            get_velocity_new_const(), get_density_new_const(), get_tracer_new_const(),
                            GetVecOfPtrs(u_mac), GetVecOfPtrs(v_mac),
                            GetVecOfPtrs(w_mac), 
                            {}, {}, new_time);

    // *************************************************************************************
    // Compute viscosity / diffusive coefficients
    // *************************************************************************************
    compute_viscosity(GetVecOfPtrs(vel_eta), GetVecOfPtrs(tra_eta),
                      get_density_new_const(), get_velocity_new_const(), get_tracer_new_const(),
                      new_time, 1);

    // Here we create divtau of the (n+1,*) state that was computed in the predictor;
    //      we use this laps only if DiffusionType::Explicit
    if (m_diff_type == DiffusionType::Explicit) {
        get_diffusion_tensor_op()->compute_divtau(get_divtau_new(),
                                                  get_velocity_new_const(),
                                                  get_density_new_const(),
                                                  GetVecOfConstPtrs(vel_eta));
        if (m_advect_tracer) {
            get_diffusion_scalar_op()->compute_laps(get_laps_new(),
                                                    get_tracer_new_const(),
                                                    get_density_new_const(),
                                                    GetVecOfConstPtrs(tra_eta));
        }
    }

    // *************************************************************************************
    // Define local variables for lambda to capture.
    // *************************************************************************************
    Real l_dt = m_dt;
    bool l_constant_density = m_constant_density;
    int l_ntrac = (m_advect_tracer) ? m_ntrac : 0;

    // *************************************************************************************
    // Update density first
    // *************************************************************************************
    if (l_constant_density)
    {
        for (int lev = 0; lev <= finest_level; lev++)
            MultiFab::Copy(density_nph[lev], m_leveldata[lev]->density_o, 0, 0, 1, 0);
    } else {
        for (int lev = 0; lev <= finest_level; lev++)
        {
            auto& ld = *m_leveldata[lev];
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(ld.velocity,TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                Box const& bx = mfi.tilebox();
                Array4<Real const> const& rho_o  = ld.density_o.const_array(mfi);
                Array4<Real> const& rho_n        = ld.density.array(mfi);
                Array4<Real> const& rho_nph      = density_nph[lev].array(mfi);
                Array4<Real const> const& drdt_o = ld.conv_density_o.const_array(mfi);
                Array4<Real const> const& drdt   = ld.conv_density.const_array(mfi);

                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                        rho_n  (i,j,k) = rho_o(i,j,k) + l_dt * 0.5*(drdt(i,j,k)+drdt_o(i,j,k));
                        rho_nph(i,j,k) = 0.5 * (rho_o(i,j,k) + rho_n(i,j,k));
                });
            } // mfi
        } // lev
    } // not constant density

    // *************************************************************************************
    // Compute the tracer forcing terms (forcing for (rho s), not for s)
    // *************************************************************************************
    if (m_advect_tracer)
        compute_tra_forces(GetVecOfPtrs(tra_forces),  GetVecOfConstPtrs(density_nph));

    // *************************************************************************************
    // Update the tracer next (note that dtdt already has rho in it)
    // (rho trac)^new = (rho trac)^old + dt * (
    //                   div(rho trac u) + div (mu grad trac) + rho * f_t
    // *************************************************************************************
    if (m_advect_tracer)
    {
        for (int lev = 0; lev <= finest_level; lev++)
        {
            auto& ld = *m_leveldata[lev];
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(ld.tracer,TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                Box const& bx = mfi.tilebox();
                Array4<Real const> const& tra_o   = ld.tracer_o.const_array(mfi);
                Array4<Real const> const& rho_o   = ld.density_o.const_array(mfi);
                Array4<Real      > const& tra     = ld.tracer.array(mfi);
                Array4<Real const> const& rho     = ld.density.const_array(mfi);
                Array4<Real const> const& rho_nph = density_nph[lev].const_array(mfi);
                Array4<Real const> const& dtdt_o  = ld.conv_tracer_o.const_array(mfi);
                Array4<Real const> const& dtdt    = ld.conv_tracer.const_array(mfi);
                Array4<Real const> const& tra_f   = (l_ntrac > 0) ? tra_forces[lev].const_array(mfi)
                                                                : Array4<Real const>{};

                if (m_diff_type == DiffusionType::Explicit) 
                {
                    Array4<Real const> const& laps_o = (l_ntrac > 0) ? ld.laps_o.const_array(mfi)
                                                                     : Array4<Real const>{};
                    Array4<Real const> const& laps   = (l_ntrac > 0) ? ld.laps.const_array(mfi)
                                                                     : Array4<Real const>{};
                    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                    {
                        for (int n = 0; n < l_ntrac; ++n) 
                        {
                            tra(i,j,k,n) = rho_o(i,j,k)*tra_o(i,j,k) + l_dt *
                                ( 0.5*(  dtdt(i,j,k,n) + dtdt_o(i,j,k,n))
                                 +0.5*(laps_o(i,j,k,n) +   laps(i,j,k,n))
                                   +    tra_f(i,j,k,n) );

                            tra(i,j,k,n) /= rho(i,j,k);
                        }
                    });
                } 
                else if (m_diff_type == DiffusionType::Crank_Nicolson) 
                {
                    Array4<Real const> const& laps_o = (l_ntrac > 0) ? ld.laps_o.const_array(mfi)
                                                                     : Array4<Real const>{};
                    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                    {
                        for (int n = 0; n < l_ntrac; ++n) 
                        {
                            tra(i,j,k,n) = rho_o(i,j,k)*tra_o(i,j,k) + l_dt *
                                ( 0.5*(  dtdt(i,j,k,n) + dtdt_o(i,j,k,n))
                                 +0.5*(laps_o(i,j,k,n)                  )
                                   +    tra_f(i,j,k,n) );

                            tra(i,j,k,n) /= rho(i,j,k);
                        }
                    });
                } 
                else if (m_diff_type == DiffusionType::Implicit) 
                {
                    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                    {
                        for (int n = 0; n < l_ntrac; ++n) 
                        {
                            tra(i,j,k,n) = rho_o(i,j,k)*tra_o(i,j,k) + l_dt *
                                ( 0.5*( dtdt(i,j,k,n)+dtdt_o(i,j,k,n))
                                   +   tra_f(i,j,k,n) );

                            tra(i,j,k,n) /= rho(i,j,k);
                        }
                    });
                }
            } // mfi
        } // lev
    } // if (m_advect_tracer)

    // *************************************************************************************
    // Solve diffusion equation for tracer
    // *************************************************************************************
    if ( m_advect_tracer &&
        (m_diff_type == DiffusionType::Crank_Nicolson || m_diff_type == DiffusionType::Implicit) )
    {
        const int ng_diffusion = 1;
        for (int lev = 0; lev <= finest_level; ++lev) 
            fillphysbc_tracer(lev, new_time, m_leveldata[lev]->tracer, ng_diffusion);

        Real dt_diff = (m_diff_type == DiffusionType::Implicit) ? m_dt : 0.5*m_dt;
        get_diffusion_scalar_op()->diffuse_scalar(get_tracer_new(),
                                                  get_density_new(),
                                                  GetVecOfConstPtrs(tra_eta),
                                                  dt_diff);
    }

    // *************************************************************************************
    // Define the forcing terms to use in the final update (using half-time density)
    // *************************************************************************************
    compute_vel_forces(GetVecOfPtrs(vel_forces), get_velocity_new_const(), 
                       GetVecOfConstPtrs(density_nph), 
                       get_tracer_old_const(), get_tracer_new_const());

    // *************************************************************************************
    // Update velocity
    // *************************************************************************************
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
            Array4<Real const> const& vel_o = ld.velocity_o.const_array(mfi);
            Array4<Real const> const& dvdt = ld.conv_velocity.const_array(mfi);
            Array4<Real const> const& dvdt_o = ld.conv_velocity_o.const_array(mfi);
            Array4<Real const> const& vel_f = vel_forces[lev].const_array(mfi);

            if (m_diff_type == DiffusionType::Implicit) 
            {
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                        vel(i,j,k,idim) = vel_o(i,j,k,idim) + l_dt*
                            (0.5*(dvdt_o(i,j,k,idim)+dvdt(i,j,k,idim))+vel_f(i,j,k,idim));
                    }
                });
            } 
            else if (m_diff_type == DiffusionType::Crank_Nicolson)
            { 
                Array4<Real const> const& divtau_o = ld.divtau_o.const_array(mfi);
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                        vel(i,j,k,idim) = vel_o(i,j,k,idim) + l_dt*
                            (0.5*(dvdt_o(i,j,k,idim)+dvdt(i,j,k,idim))+vel_f(i,j,k,idim) 
                                +divtau_o(i,j,k,idim));
                    }
                });
            }
            else if (m_diff_type == DiffusionType::Explicit)
            {
                Array4<Real const> const& divtau_o = ld.divtau_o.const_array(mfi);
                Array4<Real const> const& divtau   = ld.divtau.const_array(mfi);
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                        vel(i,j,k,idim) = vel_o(i,j,k,idim) + l_dt*
                            (0.5*(dvdt_o(i,j,k,idim)+dvdt(i,j,k,idim)
                                +divtau_o(i,j,k,idim)+divtau(i,j,k,idim))
                             +vel_f(i,j,k,idim));
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
        for (int lev = 0; lev <= finest_level; ++lev)  
            fillphysbc_velocity(lev, new_time, m_leveldata[lev]->velocity, ng_diffusion);

        Real dt_diff = (m_diff_type == DiffusionType::Implicit) ? m_dt : 0.5*m_dt;
        get_diffusion_tensor_op()->diffuse_velocity(get_velocity_new(),
                                                    get_density_new(),
                                                    GetVecOfConstPtrs(vel_eta),
                                                    dt_diff);
    }

    // **********************************************************************************************
    // 
    // Project velocity field, update pressure
    bool incremental = false;
    ApplyProjection(GetVecOfConstPtrs(density_nph),new_time, m_dt, incremental);

    // **********************************************************************************************
    // 
    // Over-write velocity in cells with vfrac < 1e-4
    // 
    // **********************************************************************************************
    incflo_correct_small_cells(get_velocity_new(),
                               GetVecOfConstPtrs(u_mac), GetVecOfConstPtrs(v_mac),
                               GetVecOfConstPtrs(w_mac));
}
