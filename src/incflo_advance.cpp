#include <incflo.H>
#include "PlaneAveraging.H"
#include "Physics.H"
#include <cmath>
#include "field_ops.H"

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
        m_t_old[lev] = m_time.current_time();
        m_t_new[lev] = m_time.new_time();
    }

    m_repo.get_field("velocity").advance_states();
    m_repo.get_field("density").advance_states();
    m_repo.get_field("tracer").advance_states();

    m_repo.get_field("velocity",amr_wind::FieldState::Old).fillpatch(m_time.current_time());
    m_repo.get_field("density",amr_wind::FieldState::Old).fillpatch(m_time.current_time());
    m_repo.get_field("tracer",amr_wind::FieldState::Old).fillpatch(m_time.current_time());

    for (auto& pp: m_physics)
        pp->pre_advance_work();
    
    ApplyPredictor();

    if (!m_use_godunov) {

        m_repo.get_field("velocity",amr_wind::FieldState::New).fillpatch(m_time.current_time());
        m_repo.get_field("density",amr_wind::FieldState::New).fillpatch(m_time.current_time());
        m_repo.get_field("tracer",amr_wind::FieldState::New).fillpatch(m_time.current_time());

        ApplyCorrector();
    }

    if (m_verbose > 1)
    {
        amrex::Print() << "End of time step: " << std::endl;
        PrintMaxValues(m_time.new_time());
    }

#if 0
    if (m_test_tracer_conservation) {
        amrex::Print() << "Sum tracer volume wgt = " << m_time.current_time()+dt << "   " << volWgtSum(0,*tracer[0],0) << std::endl;
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
    BL_PROFILE("incflo::ApplyPredictor");

    // We use the new time value for things computed on the "*" state
    Real new_time = m_time.new_time();

    if (m_verbose > 2)
    {
        amrex::Print() << "Before predictor step:" << std::endl;
        PrintMaxValues(new_time);
    }

    auto& velocity_old = m_repo.get_field("velocity", amr_wind::FieldState::Old);
    auto& velocity_new = m_repo.get_field("velocity", amr_wind::FieldState::New);
    auto& density_old = m_repo.get_field("density", amr_wind::FieldState::Old);
    auto& density_new = m_repo.get_field("density", amr_wind::FieldState::New);
    auto& tracer_old = m_repo.get_field("tracer", amr_wind::FieldState::Old);
    auto& tracer_new = m_repo.get_field("tracer", amr_wind::FieldState::New);

    // only the old states are used in predictor
    auto& divtau = m_repo.get_field("divtau", amr_wind::FieldState::Old);
    auto& laps = m_repo.get_field("laps", amr_wind::FieldState::Old);
    auto& conv_velocity = m_repo.get_field("conv_velocity", amr_wind::FieldState::Old);
    auto& conv_density = m_repo.get_field("conv_density", amr_wind::FieldState::Old);
    auto& conv_tracer = m_repo.get_field("conv_tracer", amr_wind::FieldState::Old);

    auto& u_mac = m_repo.get_field("u_mac");
    auto& v_mac = m_repo.get_field("v_mac");
    auto& w_mac = m_repo.get_field("w_mac");

    if (nghost_mac() > 0) {
        for (int lev = 0; lev <= finest_level; ++lev) {
            u_mac(lev).setBndry(0.0);
            v_mac(lev).setBndry(0.0);
            w_mac(lev).setBndry(0.0);
        }
    }

    // Allocate scratch space for half time density and tracer
    auto density_nph = m_repo.create_scratch_field(1,1);
    auto tracer_nph = m_repo.create_scratch_field(m_ntrac,1);

    // Allocate scratch space for momentum eqns forces and viscosity
    auto vel_forces = m_repo.create_scratch_field(AMREX_SPACEDIM, nghost_force());
    auto vel_eta = m_repo.create_scratch_field(1, 1);

    // Allocate scratch space for tracer eqns forces and viscosity
    // fixme should this have an if statement on it like before
    auto tra_forces = m_repo.create_scratch_field(m_ntrac, nghost_force());
    auto tra_eta = m_repo.create_scratch_field(m_ntrac, 1);

    // *************************************************************************************
    // Define the forcing terms to use in the Godunov prediction
    // *************************************************************************************
    if (m_use_godunov)
    {
        compute_vel_forces((*vel_forces).vec_ptrs(),
                           velocity_old.vec_const_ptrs(),
                           density_old.vec_const_ptrs(),
                           tracer_old.vec_const_ptrs());

        // Note this is forcing for (rho s), not for s
        if (m_advect_tracer)
           compute_tra_forces((*tra_forces).vec_ptrs(), density_old.vec_const_ptrs());
    }

    // *************************************************************************************
    // Compute viscosity / diffusive coefficients
    // *************************************************************************************
    compute_viscosity((*vel_eta).vec_ptrs(), (*tra_eta).vec_ptrs(),
                      density_old.vec_const_ptrs(), velocity_old.vec_const_ptrs(), tracer_old.vec_const_ptrs(),
                      m_time.current_time(), 1);

    // *************************************************************************************
    // Compute explicit viscous term
    // *************************************************************************************
    if (need_divtau()) {
        get_diffusion_tensor_op()->compute_divtau(get_divtau_old(),
                                                  velocity_old.vec_const_ptrs(),
                                                  density_old.vec_const_ptrs(),
                                                  (*vel_eta).vec_const_ptrs());
        if (m_use_godunov)
            amr_wind::field_ops::add(*vel_forces, divtau, 0, 0, AMREX_SPACEDIM, 0);

    }



    // *************************************************************************************
    // Compute explicit diffusive terms
    // *************************************************************************************
    if (m_advect_tracer && need_divtau()) {
        get_diffusion_scalar_op()->compute_laps(laps.vec_ptrs(),
                                                tracer_old.vec_const_ptrs(),
                                                density_old.vec_const_ptrs(),
                                                (*tra_eta).vec_const_ptrs());
        if (m_use_godunov)
            amr_wind::field_ops::add(*tra_forces, laps, 0, 0, AMREX_SPACEDIM, 0);

    }

    if (m_use_godunov) {
        // fixme
        fillpatch_force(m_time.current_time(), (*vel_forces).vec_ptrs(), nghost_force());
        if (m_advect_tracer) {
            fillpatch_force(m_time.current_time(), (*tra_forces).vec_ptrs(), nghost_force());
        }
    }

    // *************************************************************************************
    // if ( m_use_godunov) Compute the explicit advective terms R_u^(n+1/2), R_s^(n+1/2) and R_t^(n+1/2)
    // if (!m_use_godunov) Compute the explicit advective terms R_u^n      , R_s^n       and R_t^n
    // Note that "get_conv_tracer_old" returns div(rho u tracer)
    // *************************************************************************************
    compute_convective_term(conv_velocity.vec_ptrs(),
                            conv_density.vec_ptrs(),
                            conv_tracer.vec_ptrs(),
                            velocity_old.vec_const_ptrs(),
                            density_old.vec_const_ptrs(),
                            tracer_old.vec_const_ptrs(),
                            u_mac.vec_ptrs(),
                            v_mac.vec_ptrs(),
                            w_mac.vec_ptrs(),
                            (*vel_forces).vec_const_ptrs(),
                            (*tra_forces).vec_const_ptrs(),
                            m_time.current_time());

    // *************************************************************************************
    // Define local variables for lambda to capture.
    // *************************************************************************************
    Real l_dt = m_time.deltaT();
    bool l_constant_density = m_constant_density;

    // *************************************************************************************
    // Update density first
    // *************************************************************************************
    if (l_constant_density)
    {
        amr_wind::field_ops::copy(*density_nph, density_old, 0, 0, 1, 1);
    }
    else
    {

        for (int lev = 0; lev <= finest_level; lev++)
        {

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(density_old(lev),TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                Box const& bx = mfi.tilebox();
                Array4<Real  const> const& rho_o  = density_old(lev).const_array(mfi);
                Array4<Real> const& rho_new       = density_new(lev).array(mfi);
                Array4<Real> const& rho_nph       = (*density_nph)(lev).array(mfi);
                Array4<Real const> const& drdt    = conv_density(lev).const_array(mfi);

                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    rho_new(i,j,k) = rho_o(i,j,k) + l_dt * drdt(i,j,k);
                    rho_nph(i,j,k) = 0.5 * (rho_o(i,j,k) + rho_new(i,j,k));
                });
            } // mfi
        } // lev

    } // not constant density

    // *************************************************************************************
    // Compute (or if Godunov, re-compute) the tracer forcing terms (forcing for (rho s), not for s)
    // *************************************************************************************
    if (m_advect_tracer)
       compute_tra_forces((*tra_forces).vec_ptrs(), (*density_nph).vec_const_ptrs());

    // *************************************************************************************
    // Update the tracer next
    // *************************************************************************************
    int l_ntrac = (m_advect_tracer) ? m_ntrac : 0;

    if (m_advect_tracer)
    {
        for (int lev = 0; lev <= finest_level; lev++)
        {

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(tracer_old(lev),TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                Box const& bx = mfi.tilebox();
                Array4<Real const> const& tra_o   = tracer_old(lev).const_array(mfi);
                Array4<Real const> const& rho_o   = density_old(lev).const_array(mfi);
                Array4<Real> const& tra           = tracer_new(lev).array(mfi);
                Array4<Real const> const& rho     = density_new(lev).const_array(mfi);
                Array4<Real const> const& dtdt_o  = conv_tracer(lev).const_array(mfi);
                Array4<Real const> const& tra_f   = (*tra_forces)(lev).const_array(mfi);

                if (m_diff_type == DiffusionType::Explicit)
                {
                    Array4<Real const> const& laps_o = laps(lev).const_array(mfi);
                    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                    {
                        // (rho trac)^new = (rho trac)^old + dt * (
                        //                   div(rho trac u) + div (mu grad trac) + rho * f_t
                        for (int n = 0; n < l_ntrac; ++n)
                        {
                            tra(i,j,k,n) = rho_o(i,j,k)*tra_o(i,j,k,n) + l_dt *
                                ( dtdt_o(i,j,k,n) + tra_f(i,j,k,n) + laps_o(i,j,k,n) );

                            tra(i,j,k,n) /= rho(i,j,k);
                        }
                    });
                }
                else if (m_diff_type == DiffusionType::Crank_Nicolson)
                {
                    Array4<Real const> const& laps_o = laps(lev).const_array(mfi);
                    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                    {
                        for (int n = 0; n < l_ntrac; ++n)
                        {
                            tra(i,j,k,n) = rho_o(i,j,k)*tra_o(i,j,k,n) + l_dt *
                                ( dtdt_o(i,j,k,n) + tra_f(i,j,k,n) + 0.5 * laps_o(i,j,k,n) );

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
                            tra(i,j,k,n) = rho_o(i,j,k)*tra_o(i,j,k,n) + l_dt *
                                ( dtdt_o(i,j,k,n) + tra_f(i,j,k,n) );

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
        IntVect ng_diffusion(1);
        tracer_new.fillphysbc(new_time, ng_diffusion);

        Real dt_diff = (m_diff_type == DiffusionType::Implicit) ? m_time.deltaT() : 0.5*m_time.deltaT();
        get_diffusion_scalar_op()->diffuse_scalar(tracer_new.vec_ptrs(),
                                                  density_new.vec_ptrs(),
                                                  (*tra_eta).vec_const_ptrs(),
                                                  dt_diff);
    } // if (m_advect_tracer)

    
    // *************************************************************************************
    // Update tracer at n+1/2
    // *************************************************************************************
    if (m_advect_tracer)
        amr_wind::field_ops::lincomb(*tracer_nph, 0.5, tracer_old, 0, 0.5, tracer_new, 0, 0, m_ntrac, 1);
    else
        amr_wind::field_ops::copy(*tracer_nph, tracer_old, 0, 0, m_ntrac, 1);

    
    
    // *************************************************************************************
    // Define (or if use_godunov, re-define) the forcing terms, without the viscous terms
    //    and using the half-time density
    // *************************************************************************************
    compute_vel_forces((*vel_forces).vec_ptrs(),
                       velocity_old.vec_const_ptrs(),
                       (*density_nph).vec_const_ptrs(),
                       (*tracer_nph).vec_const_ptrs());
    
    // *************************************************************************************
    // Update the velocity
    // *************************************************************************************
    for (int lev = 0; lev <= finest_level; lev++)
    {

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(velocity_new(lev),TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            Box const& bx = mfi.tilebox();
            Array4<Real> const& vel = velocity_new(lev).array(mfi);
            Array4<Real const> const& dvdt = conv_velocity(lev).const_array(mfi);
            Array4<Real const> const& vel_f = (*vel_forces)(lev).const_array(mfi);

            if (m_diff_type == DiffusionType::Implicit) {

                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    vel(i,j,k,0) += l_dt*(dvdt(i,j,k,0)+vel_f(i,j,k,0));
                    vel(i,j,k,1) += l_dt*(dvdt(i,j,k,1)+vel_f(i,j,k,1));
                    vel(i,j,k,2) += l_dt*(dvdt(i,j,k,2)+vel_f(i,j,k,2));
                });
            }
            else if (m_diff_type == DiffusionType::Crank_Nicolson)
            {

                Array4<Real const> const& divtau_o = divtau(lev).const_array(mfi);
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    vel(i,j,k,0) += l_dt*(dvdt(i,j,k,0)+vel_f(i,j,k,0)+0.5*divtau_o(i,j,k,0));
                    vel(i,j,k,1) += l_dt*(dvdt(i,j,k,1)+vel_f(i,j,k,1)+0.5*divtau_o(i,j,k,1));
                    vel(i,j,k,2) += l_dt*(dvdt(i,j,k,2)+vel_f(i,j,k,2)+0.5*divtau_o(i,j,k,2));
                });
            }
            else if (m_diff_type == DiffusionType::Explicit)
            {
                Array4<Real const> const& divtau_o = divtau(lev).const_array(mfi);
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    vel(i,j,k,0) += l_dt*(dvdt(i,j,k,0)+vel_f(i,j,k,0)+divtau_o(i,j,k,0));
                    vel(i,j,k,1) += l_dt*(dvdt(i,j,k,1)+vel_f(i,j,k,1)+divtau_o(i,j,k,1));
                    vel(i,j,k,2) += l_dt*(dvdt(i,j,k,2)+vel_f(i,j,k,2)+divtau_o(i,j,k,2));
                });
            }
        } // mfi
    } // lev

    // *************************************************************************************
    // Solve diffusion equation for u* but using eta_old at old time
    // *************************************************************************************
    if (m_diff_type == DiffusionType::Crank_Nicolson || m_diff_type == DiffusionType::Implicit)
    {
        IntVect ng_diffusion(1);
        velocity_new.fillphysbc(new_time, ng_diffusion);
        density_new.fillphysbc(new_time, ng_diffusion);

        Real dt_diff = (m_diff_type == DiffusionType::Implicit) ? m_time.deltaT() : 0.5*m_time.deltaT();
        get_diffusion_tensor_op()->diffuse_velocity(velocity_new.vec_ptrs(),
                                                    density_new.vec_const_ptrs(),
                                                    (*vel_eta).vec_const_ptrs(),
                                                    dt_diff);
    }

    // **********************************************************************************************
    //
    // Project velocity field, update pressure
    //
    // **********************************************************************************************
    ApplyProjection((*density_nph).vec_const_ptrs(), new_time, m_time.deltaT(), incremental_projection);

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
    Real new_time = m_time.new_time();

    if (m_verbose > 2)
    {
        amrex::Print() << "Before corrector step:" << std::endl;
        PrintMaxValues(new_time);
    }

    // *************************************************************************************
    // Allocate space for half-time density
    // *************************************************************************************
    Vector<MultiFab> density_nph;
    Vector<MultiFab> tracer_nph;
    for (int lev = 0; lev <= finest_level; ++lev) {
        density_nph.emplace_back(grids[lev], dmap[lev], 1, 0, MFInfo(), Factory(lev));
        if (m_ntrac) tracer_nph.emplace_back(grids[lev], dmap[lev], m_ntrac, 1, MFInfo(), Factory(lev));

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
                            m_repo.get_field("u_mac").vec_ptrs(),
                            m_repo.get_field("v_mac").vec_ptrs(),
                            m_repo.get_field("w_mac").vec_ptrs(),
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
    Real l_dt = m_time.deltaT();
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
                            tra(i,j,k,n) = rho_o(i,j,k)*tra_o(i,j,k,n) + l_dt *
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
                            tra(i,j,k,n) = rho_o(i,j,k)*tra_o(i,j,k,n) + l_dt *
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
                            tra(i,j,k,n) = rho_o(i,j,k)*tra_o(i,j,k,n) + l_dt *
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

        Real dt_diff = (m_diff_type == DiffusionType::Implicit) ? m_time.deltaT() : 0.5*m_time.deltaT();
        get_diffusion_scalar_op()->diffuse_scalar(get_tracer_new(),
                                                  get_density_new(),
                                                  GetVecOfConstPtrs(tra_eta),
                                                  dt_diff);
    }

    // *************************************************************************************
    // Update tracer at n+1/2
    // *************************************************************************************
    if (!m_advect_tracer)
    {
        for (int lev = 0; lev <= finest_level; lev++){
            MultiFab::Copy(tracer_nph[lev], m_leveldata[lev]->tracer_o, 0, 0, m_ntrac, 1);
        }
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
                Array4<Real const> const& tracer_o  = ld.tracer_o.const_array(mfi);
                Array4<Real const> const& tracer_n  = ld.tracer.const_array(mfi);
                Array4<Real> const& tra_nph   = tracer_nph[lev].array(mfi);

                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    tra_nph(i,j,k) = 0.5 * (tracer_o(i,j,k) + tracer_n(i,j,k));
                });
            } // mfi
        } // lev

    }
    
    // *************************************************************************************
    // Define the forcing terms to use in the final update (using half-time density)
    // *************************************************************************************
    compute_vel_forces(GetVecOfPtrs(vel_forces), get_velocity_new_const(),
                       GetVecOfConstPtrs(density_nph),
                       GetVecOfConstPtrs(tracer_nph));

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
        {
            fillphysbc_velocity(lev, new_time, m_leveldata[lev]->velocity, ng_diffusion);
            fillphysbc_density (lev, new_time, m_leveldata[lev]->density , ng_diffusion);
        }
        Real dt_diff = (m_diff_type == DiffusionType::Implicit) ? m_time.deltaT() : 0.5*m_time.deltaT();
        get_diffusion_tensor_op()->diffuse_velocity(get_velocity_new(),
                                                    get_density_new_const(),
                                                    GetVecOfConstPtrs(vel_eta),
                                                    dt_diff);
    }

    // **********************************************************************************************
    //
    // Project velocity field, update pressure
    bool incremental = false;
    ApplyProjection(GetVecOfConstPtrs(density_nph),new_time, m_time.deltaT(), incremental);

}
