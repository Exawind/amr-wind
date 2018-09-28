#include <AMReX_Array.H>
#include <AMReX_BC_TYPES.H>
#include <AMReX_BLassert.H>
#include <AMReX_Box.H>
#include <AMReX_EBMultiFabUtil.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_VisMF.H>
#include <incflo_icbc_F.H>
#include <incflo_level.H>
#include <incflo_mac_F.H>
#include <incflo_proj_F.H>

void incflo_level::Advance(
	int lev, int nstep, int steady_state, Real& dt, Real& prev_dt, Real time, Real stop_time)
{
	AMREX_ALWAYS_ASSERT(lev == 0);

	BL_PROFILE_REGION_START("incflo::Advance");
	BL_PROFILE("incflo::Advance");

	amrex::Print() << "\n ============   NEW TIME STEP   ============ \n";

	// Extrapolate boundary values for density and volume fraction
	fill_mf_bc(lev, *mu[lev]);

	// Fill ghost nodes and reimpose boundary conditions
	incflo_set_scalar_bcs(lev);
	incflo_set_velocity_bcs(lev, 0);

	// Start loop: if we are not seeking a steady state solution,
	// the loop will execute only once
	int keep_looping = 1;
	int iter = 1;
	do
	{
        // Compute time step size
		incflo_compute_dt(lev, time, stop_time, steady_state, dt);

		if(steady_state)
		{
			amrex::Print() << "\n   Iteration " << iter << " with dt = " << dt << "\n" << std::endl;
		}
		else
		{
			amrex::Print() << "\n   Step " << nstep + 1 << ": from old_time " << time
						   << " to new time " << time + dt << " with dt = " << dt << "\n"
						   << std::endl;
		}

		// Backup field variable to old
		MultiFab::Copy(*p_o[lev], *p[lev], 0, 0, p[lev]->nComp(), p_o[lev]->nGrow());
		MultiFab::Copy(*ro_o[lev], *ro[lev], 0, 0, ro[lev]->nComp(), ro_o[lev]->nGrow());
		MultiFab::Copy(*vel_o[lev], *vel[lev], 0, 0, vel[lev]->nComp(), vel_o[lev]->nGrow());

		// Time integration step
		//
		// Create temporary multifabs to hold the old-time conv and divtau
		//    so we don't have to re-compute them in the corrector
		MultiFab conv_old(grids[lev], dmap[lev], 3, 0, MFInfo(), *ebfactory[lev]);
		MultiFab divtau_old(grids[lev], dmap[lev], 3, 0, MFInfo(), *ebfactory[lev]);

		// Predictor step
		bool proj_2 = true;
		incflo_apply_predictor(lev, conv_old, divtau_old, dt, proj_2);

		// Print info about predictor step
        if(verbose > 0)
		{
			amrex::Print() << "\nAfter predictor step:\n";
			incflo_print_max_vel(lev);
			incflo_compute_diveu(lev);
			amrex::Print() << "max(abs(diveu)) = " << incflo_norm0(diveu, lev, 0) << "\n";
		}

		// Corrector step
		proj_2 = true;
		incflo_apply_corrector(lev, conv_old, divtau_old, dt, proj_2);

		// Print info about corrector step
        if(verbose > 0)
		{
			amrex::Print() << "\nAfter corrector step:\n";
			incflo_print_max_vel(lev);
			incflo_compute_diveu(lev);
			amrex::Print() << "max(abs(diveu)) = " << incflo_norm0(diveu, lev, 0) << "\n";
		}

		//
		// Check whether to exit the loop or not
		//
		if(steady_state)
			keep_looping = !steady_state_reached(lev, dt);
		else
			keep_looping = 0;

		// Update interations count
		++iter;
	} while(keep_looping);

    prev_dt = dt;

	BL_PROFILE_REGION_STOP("incflo::Advance");
}

void incflo_level::incflo_compute_dt(int lev, Real time, Real stop_time, int steady_state, Real& dt)
{
	// DT is always computed even for fixed dt, so we can
	// issue a warning if fixed dt does not satisfy CFL condition.
	Real dt_new = dt;

	// Compute dt for this time step
	Real umax = incflo_norm0(vel, lev, 0);
	Real vmax = incflo_norm0(vel, lev, 1);
	Real wmax = incflo_norm0(vel, lev, 2);
	Real romax = incflo_norm0(ro, lev, 0);
	Real mumax = incflo_norm0(mu, lev, 0);

	Real gradp0max[3];

	if(nodal_pressure == 1)
	{
		for(MFIter mfi(*vel[lev], true); mfi.isValid(); ++mfi)
		{
			// Cell-centered tilebox
			Box bx = mfi.tilebox();

			compute_gradp0_max(bx.loVect(),
							   bx.hiVect(),
							   BL_TO_FORTRAN_ANYD((*p0[lev])[mfi]),
							   gradp0max,
							   geom[lev].CellSize(),
							   &nodal_pressure);
		}
	}
	else
	{
		gradp0max[0] = incflo_norm0(gp0, lev, 0);
		gradp0max[1] = incflo_norm0(gp0, lev, 1);
		gradp0max[2] = incflo_norm0(gp0, lev, 2);
	}

	ParallelDescriptor::ReduceRealMax(gradp0max[0]);
	ParallelDescriptor::ReduceRealMax(gradp0max[1]);
	ParallelDescriptor::ReduceRealMax(gradp0max[2]);

	compute_new_dt(&umax,
				   &vmax,
				   &wmax,
				   &romax,
				   &mumax,
				   gradp0max,
				   geom[lev].CellSize(),
				   &cfl,
				   &steady_state,
				   &time,
				   &stop_time,
				   &dt_new);

	if(fixed_dt > 0.)
	{
		if(dt_new < fixed_dt)
		{
			amrex::Print() << "WARNING: fixed_dt does not satisfy CFL condition: "
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

void incflo_level::incflo_project_velocity(int lev)
{
	// Project velocity field to make sure initial velocity is divergence-free

	amrex::Print() << "Initial projection:\n";

	// Need to add this call here so that the MACProjection internal arrays
	//  are allocated so that the cell-centered projection can use the MAC
	//  data structures and set_velocity_bcs routine
	mac_projection->update_internals();

	bool proj_2 = true;
	Real dummy_dt = 1.0;
	incflo_apply_projection(lev, dummy_dt, proj_2);

	// We initialize p and gp back to zero (p0 may still be still non-zero)
	p[lev]->setVal(0.0);
	 gp[lev]->setVal(0.0);
}

void incflo_level::incflo_initial_iterations(int lev, Real dt, Real stop_time, int steady_state)
{
	// Fill ghost cells
	incflo_set_scalar_bcs(lev);
	incflo_set_velocity_bcs(lev, 0);

	// Copy vel into vel_o
	MultiFab::Copy(*vel_o[lev], *vel[lev], 0, 0, vel[lev]->nComp(), vel_o[lev]->nGrow());

	Real time = 0.0;
	incflo_compute_dt(lev, time, stop_time, steady_state, dt);

	amrex::Print() << "Doing initial pressure iterations with dt = " << dt << std::endl;

	//  Create temporary multifabs to hold conv and divtau
	MultiFab conv(grids[lev], dmap[lev], 3, 0, MFInfo(), *ebfactory[lev]);
	MultiFab divtau(grids[lev], dmap[lev], 3, 0, MFInfo(), *ebfactory[lev]);

	for(int iter = 0; iter < 3; ++iter)
	{
		amrex::Print() << "\n In initial_iterations: iter = " << iter << "\n";

		bool proj_2 = false;
		incflo_apply_predictor(lev, conv, divtau, dt, proj_2);

		// Replace vel by the original values
		MultiFab::Copy(*vel[lev], *vel_o[lev], 0, 0, vel[lev]->nComp(), vel[lev]->nGrow());
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
void incflo_level::incflo_apply_predictor(
	int lev, MultiFab& conv_old, MultiFab& divtau_old, amrex::Real dt, bool proj_2)
{
	// Compute the explicit advective term R_u^n
	incflo_compute_ugradu_predictor(lev, conv_old, vel_o);

	// If explicit_diffusion == true  then we compute the full diffusive terms
	// here
	// If explicit_diffusion == false then we compute only the off-diagonal terms
	// here
	incflo_compute_divtau(lev, divtau_old, vel_o);

	// First add the convective term
	MultiFab::Saxpy(*vel[lev], dt, conv_old, 0, 0, 3, 0);

	// Add the diffusion terms (either all if explicit_diffusion == true or just
	// the
	//    off-diagonal terms if explicit_diffusion == false)
	MultiFab::Saxpy(*vel[lev], dt, divtau_old, 0, 0, 3, 0);

	// Add the forcing terms
	incflo_apply_forcing_terms(lev, dt, vel);

	// Convert velocities to momenta
	for(int n = 0; n < 3; n++)
		MultiFab::Multiply(*vel[lev], (*ro[lev]), 0, n, 1, vel[lev]->nGrow());

	// Add (-dt grad p to momenta)
	MultiFab::Saxpy(*vel[lev], -dt, *gp[lev], 0, 0, 3, vel[lev]->nGrow());
	MultiFab::Saxpy(*vel[lev], -dt, *gp0[lev], 0, 0, 3, vel[lev]->nGrow());

	// Convert momenta back to velocities
	for(int n = 0; n < 3; n++)
		MultiFab::Divide(*vel[lev], (*ro[lev]), 0, n, 1, vel[lev]->nGrow());

	// If doing implicit diffusion, solve here for u^*
	if(!explicit_diffusion)
		incflo_diffuse_velocity(lev, dt);

	// Project velocity field
	incflo_apply_projection(lev, dt, proj_2);
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
void incflo_level::incflo_apply_corrector(
	int lev, MultiFab& conv_old, MultiFab& divtau_old, amrex::Real dt, bool proj_2)
{
	BL_PROFILE("incflo_level::incflo_apply_corrector");

	MultiFab conv(grids[lev], dmap[lev], 3, 0);
	MultiFab divtau(grids[lev], dmap[lev], 3, 0, MFInfo(), *ebfactory[lev]);

	// Compute the explicit advective term R_u^*
	incflo_compute_ugradu_corrector(lev, conv, vel);

	// If explicit_diffusion == true  then we compute the full diffusive terms
	// here
	// If explicit_diffusion == false then we compute only the off-diagonal terms
	// here
	incflo_compute_divtau(lev, divtau, vel);

	// Define u = u_o + dt/2 (R_u^* + R_u^n)
	MultiFab::LinComb(*vel[lev], 1.0, *vel_o[lev], 0, dt / 2.0, conv, 0, 0, 3, 0);
	MultiFab::Saxpy(*vel[lev], dt / 2.0, conv_old, 0, 0, 3, 0);

	// Add the diffusion terms (either all if explicit_diffusion == true or just
	// the
	//    off-diagonal terms if explicit_diffusion == false)
	MultiFab::Saxpy(*vel[lev], dt / 2.0, divtau, 0, 0, 3, 0);
	MultiFab::Saxpy(*vel[lev], dt / 2.0, divtau_old, 0, 0, 3, 0);

	// Add forcing terms
	incflo_apply_forcing_terms(lev, dt, vel);

	// Convert velocities to momenta
	for(int n = 0; n < 3; n++)
		MultiFab::Multiply(*vel[lev], (*ro[lev]), 0, n, 1, vel[lev]->nGrow());

	// Add (-dt grad p to momenta)
	MultiFab::Saxpy(*vel[lev], -dt, *gp[lev], 0, 0, 3, vel[lev]->nGrow());
	MultiFab::Saxpy(*vel[lev], -dt, *gp0[lev], 0, 0, 3, vel[lev]->nGrow());

	// Convert momenta back to velocities
	for(int n = 0; n < 3; n++)
		MultiFab::Divide(*vel[lev], (*ro[lev]), 0, n, 1, vel[lev]->nGrow());

	// If doing implicit diffusion, solve here for u^*
	if(!explicit_diffusion)
		incflo_diffuse_velocity(lev, dt);

	// Apply projection
	incflo_apply_projection(lev, dt, proj_2);
}

void incflo_level::incflo_add_grad_phi(int lev, amrex::Real coeff, MultiFab& this_phi)
{
	BL_PROFILE("incflo_level::incflo_add_grad_phi");

	if(nodal_pressure == 1)
	{
#ifdef _OPENMP
#pragma omp parallel
#endif
		for(MFIter mfi(*vel[lev], true); mfi.isValid(); ++mfi)
		{
			// Tilebox
			Box bx = mfi.tilebox();

			add_grad_phind(BL_TO_FORTRAN_BOX(bx),
						   BL_TO_FORTRAN_ANYD((*vel[lev])[mfi]),
						   BL_TO_FORTRAN_ANYD((*ro[lev])[mfi]),
						   BL_TO_FORTRAN_ANYD(this_phi[mfi]),
						   geom[lev].CellSize(),
						   &coeff);
		}
	}
	else
	{
#ifdef _OPENMP
#pragma omp parallel
#endif
		for(MFIter mfi(*vel[lev], true); mfi.isValid(); ++mfi)
		{
			// Tilebox
			Box bx = mfi.tilebox();

			add_grad_phicc(BL_TO_FORTRAN_BOX(bx),
						   BL_TO_FORTRAN_ANYD((*vel[lev])[mfi]),
						   BL_TO_FORTRAN_ANYD((*ro[lev])[mfi]),
						   this_phi[mfi].dataPtr(),
						   geom[lev].CellSize(),
						   &coeff);
		}
	}
}

void incflo_level::incflo_apply_forcing_terms(int lev,
											  amrex::Real dt,
											  Vector<std::unique_ptr<MultiFab>>& vel)

{
	BL_PROFILE("incflo_level::incflo_apply_forcing_terms");

	Box domain(geom[lev].Domain());

#ifdef _OPENMP
#pragma omp parallel
#endif
	for(MFIter mfi(*vel[lev], true); mfi.isValid(); ++mfi)
	{
		// Tilebox
		Box bx = mfi.tilebox();

		add_forcing(BL_TO_FORTRAN_BOX(bx),
					BL_TO_FORTRAN_ANYD((*vel[lev])[mfi]),
					BL_TO_FORTRAN_ANYD((*ro[lev])[mfi]),
					domain.loVect(),
					domain.hiVect(),
					geom[lev].CellSize(),
					&dt);
	}
}

//
// Compute div(u)
//
void incflo_level::incflo_compute_diveu(int lev)
{
	Box domain(geom[lev].Domain());

	if(nodal_pressure == 1)
	{

		// Create a temporary multifab to hold (vel)
		MultiFab vec(vel[lev]->boxArray(),
					 vel[lev]->DistributionMap(),
					 vel[lev]->nComp(),
					 vel[lev]->nGrow());

		// Fill it with (vel)
		vec.copy(*vel[lev], 0, 0, vel[lev]->nComp(), vel[lev]->nGrow(), vel[lev]->nGrow());

#ifdef _OPENMP
#pragma omp parallel
#endif
		// Extrapolate Dirichlet values to ghost cells -- but do it differently in
		// that
		//  no-slip walls are treated exactly like slip walls -- this is only
		//  relevant
		//  when going into the projection
		for(MFIter mfi(vec, true); mfi.isValid(); ++mfi)
		{
			set_vec_bcs(BL_TO_FORTRAN_ANYD(vec[mfi]),
						bc_ilo.dataPtr(),
						bc_ihi.dataPtr(),
						bc_jlo.dataPtr(),
						bc_jhi.dataPtr(),
						bc_klo.dataPtr(),
						bc_khi.dataPtr(),
						domain.loVect(),
						domain.hiVect(),
						&nghost);
		}

		vec.FillBoundary(geom[lev].periodicity());

#ifdef _OPENMP
#pragma omp parallel
#endif
		for(MFIter mfi(*diveu[lev], true); mfi.isValid(); ++mfi)
		{
			// Note: this box is nodal!
			const Box& bx = mfi.tilebox();

			compute_diveund(BL_TO_FORTRAN_BOX(bx),
							BL_TO_FORTRAN_ANYD((*diveu[lev])[mfi]),
							BL_TO_FORTRAN_ANYD(vec[mfi]),
							geom[lev].CellSize());
		}
	}
	else
	{
        int extrap_dir_bcs = 1;
		incflo_set_velocity_bcs(lev, extrap_dir_bcs);
		vel[lev]->FillBoundary(geom[lev].periodicity());

        // Create face centered multifabs for and vel
        Array<std::unique_ptr<MultiFab>,AMREX_SPACEDIM> vel_fc;
        incflo_average_cc_to_fc( lev, *vel[lev], vel_fc );

        // This does not need to have correct ghost values in place
        EB_computeDivergence( *diveu[lev], GetArrOfConstPtrs(vel_fc), geom[lev] );
	}

	// Restore velocities to carry Dirichlet values on faces
	int extrap_dir_bcs = 0;
	incflo_set_velocity_bcs(lev, extrap_dir_bcs);
}

//
// Check if steady state has been reached by verifying that
//
//      max(abs( u^(n+1) - u^(n) )) < tol * dt
//      max(abs( v^(n+1) - v^(n) )) < tol * dt
//      max(abs( w^(n+1) - w^(n) )) < tol * dt
//

int incflo_level::steady_state_reached(int lev, Real dt)
{
	//
	// Count number of access
	//
	static int naccess = 0;

	//
	// Make sure velocity is up to date
	//
	incflo_set_velocity_bcs(lev, 0);

	//
	// Use temporaries to store the difference
	// between current and previous solution
	//
	MultiFab temp_vel(vel[lev]->boxArray(), vel[lev]->DistributionMap(), 3, 0);
	MultiFab::LinComb(temp_vel, 1.0, *vel[lev], 0, -1.0, *vel_o[lev], 0, 0, 3, 0);

	MultiFab tmp;

	if(nodal_pressure)
	{
		const BoxArray& nd_grid = amrex::convert(grids[lev], IntVect{1, 1, 1});
		tmp.define(nd_grid, dmap[lev], 1, 0);
	}
	else
	{
		tmp.define(grids[lev], dmap[lev], 1, 0);
	}

	MultiFab::LinComb(tmp, 1.0, *p[lev], 0, -1.0, *p_o[lev], 0, 0, 1, 0);

	Real delta_u = incflo_norm0(temp_vel, lev, 0);
	Real delta_v = incflo_norm0(temp_vel, lev, 1);
	Real delta_w = incflo_norm0(temp_vel, lev, 2);
	Real delta_p = incflo_norm0(tmp, lev, 0);

	Real tol = steady_state_tol;

	int condition1 = (delta_u < tol * dt) && (delta_v < tol * dt) && (delta_w < tol * dt);

	//
	// Second stop condition
	//
	Real du_n1 = incflo_norm1(temp_vel, lev, 0);
	Real dv_n1 = incflo_norm1(temp_vel, lev, 1);
	Real dw_n1 = incflo_norm1(temp_vel, lev, 2);
	Real dp_n1 = incflo_norm1(tmp, lev, 0);
	Real uo_n1 = incflo_norm1(vel_o, lev, 0);
	Real vo_n1 = incflo_norm1(vel_o, lev, 1);
	Real wo_n1 = incflo_norm1(vel_o, lev, 2);
	Real po_n1 = incflo_norm1(p_o, lev, 0);

	Real tmp1, tmp2, tmp3, tmp4;

	Real local_tol = 1.0e-8;

	if(uo_n1 < local_tol)
		tmp1 = 0.0;
	else
		tmp1 = du_n1 / uo_n1;

	if(vo_n1 < local_tol)
		tmp2 = 0.0;
	else
		tmp2 = dv_n1 / vo_n1;

	if(wo_n1 < local_tol)
		tmp3 = 0.0;
	else
		tmp3 = dw_n1 / wo_n1;

	if(po_n1 < local_tol)
		tmp4 = 0.0;
	else
		tmp4 = dp_n1 / po_n1;

	int condition2 = (tmp1 < tol) && (tmp2 < tol) && (tmp3 < tol); // && (tmp4 < tol);

	//
	// Print out info on steady state checks
	//
	amrex::Print() << "\nSteady state check:\n";
	amrex::Print() << "||u-uo||/||uo|| , du/dt  = " << tmp1 << " , " << delta_u / dt << "\n";
	amrex::Print() << "||v-vo||/||vo|| , dv/dt  = " << tmp2 << " , " << delta_v / dt << "\n";
	amrex::Print() << "||w-wo||/||wo|| , dw/dt  = " << tmp3 << " , " << delta_w / dt << "\n";
	amrex::Print() << "||p-po||/||po|| , dp/dt  = " << tmp4 << " , " << delta_p / dt << "\n";

	// Count # access
	naccess++;

	//
	//  Always return negative to first access. This way
	//  initial zero velocity field do not test for false positive
	//
	if(naccess == 1)
		return 0;
	else
		return condition1 || condition2;
}

//
// Set the BCs for all the variables EXCEPT pressure or velocity.
//
void incflo_level::incflo_set_scalar_bcs(int lev)
{
	BL_PROFILE("incflo_level::incflo_set_scalar_bcs()");

	Box domain(geom[lev].Domain());

#ifdef _OPENMP
#pragma omp parallel
#endif
	for(MFIter mfi(*ro[lev], true); mfi.isValid(); ++mfi)
	{
		set_scalar_bcs((*ro[lev])[mfi].dataPtr(),
					   (*mu[lev])[mfi].dataPtr(),
					   BL_TO_FORTRAN_ANYD((*lambda[lev])[mfi]),
					   bc_ilo.dataPtr(),
					   bc_ihi.dataPtr(),
					   bc_jlo.dataPtr(),
					   bc_jhi.dataPtr(),
					   bc_klo.dataPtr(),
					   bc_khi.dataPtr(),
					   domain.loVect(),
					   domain.hiVect(),
					   &nghost);
	}
	ro[lev]->FillBoundary(geom[lev].periodicity());
	mu[lev]->FillBoundary(geom[lev].periodicity());
	lambda[lev]->FillBoundary(geom[lev].periodicity());
}

//
// Set the BCs for velocity only
//
void incflo_level::incflo_set_velocity_bcs(int lev, int extrap_dir_bcs)
{
	BL_PROFILE("incflo_level::incflo_set_velocity_bcs()");

	vel[lev]->FillBoundary(geom[lev].periodicity());

	Box domain(geom[lev].Domain());

#ifdef _OPENMP
#pragma omp parallel
#endif
	for(MFIter mfi(*vel[lev], true); mfi.isValid(); ++mfi)
	{
		set_velocity_bcs(BL_TO_FORTRAN_ANYD((*vel[lev])[mfi]),
						 bc_ilo.dataPtr(),
						 bc_ihi.dataPtr(),
						 bc_jlo.dataPtr(),
						 bc_jhi.dataPtr(),
						 bc_klo.dataPtr(),
						 bc_khi.dataPtr(),
						 domain.loVect(),
						 domain.hiVect(),
						 &nghost,
						 &extrap_dir_bcs);
	}
}

//
// Fills ghost cell values of pressure appropriately for the BC type
//
void incflo_level::incflo_extrap_pressure(int lev, std::unique_ptr<amrex::MultiFab>& p)
{
	BL_PROFILE("incflo_level::incflo_extrap_pressure()");
	if(nodal_pressure == 1)
		return;

	Box domain(geom[lev].Domain());

#ifdef _OPENMP
#pragma omp parallel
#endif
	for(MFIter mfi(*p, true); mfi.isValid(); ++mfi)
	{

		extrap_pressure_to_ghost_cells(BL_TO_FORTRAN_ANYD((*p)[mfi]),
									   bc_ilo.dataPtr(),
									   bc_ihi.dataPtr(),
									   bc_jlo.dataPtr(),
									   bc_jhi.dataPtr(),
									   bc_klo.dataPtr(),
									   bc_khi.dataPtr(),
									   domain.loVect(),
									   domain.hiVect(),
									   &nghost);
	}
}

void incflo_level::check_for_nans(int lev)
{
	bool ug_has_nans = vel[lev]->contains_nan(0);
	bool vg_has_nans = vel[lev]->contains_nan(1);
	bool wg_has_nans = vel[lev]->contains_nan(2);
	bool pg_has_nans = p[lev]->contains_nan(0);

	if(ug_has_nans)
		amrex::Print() << "WARNING: u contains NaNs!!!";

	if(vg_has_nans)
		amrex::Print() << "WARNING: v contains NaNs!!!";

	if(wg_has_nans)
		amrex::Print() << "WARNING: w contains NaNs!!!";

	if(pg_has_nans)
		amrex::Print() << "WARNING: p contains NaNs!!!";
}

//
// Print the maximum values of the velocity components
//
void incflo_level::incflo_print_max_vel(int lev)
{
	amrex::Print() << "max(abs(u/v/w/p))  = " 
                   << incflo_norm0(vel, lev, 0) << "  "
				   << incflo_norm0(vel, lev, 1) << "  " 
                   << incflo_norm0(vel, lev, 2) << "  "
				   << incflo_norm0(p, lev, 0) << "  " << std::endl;
}

//
// This subroutines averages component by component
// The assumption is that cc is multicomponent
// 
void
incflo_level::incflo_average_cc_to_fc(int lev, const MultiFab& cc,
                                      Array<std::unique_ptr<MultiFab>,AMREX_SPACEDIM>& fc )
{
    AMREX_ASSERT(cc.nComp()==AMREX_SPACEDIM);
    AMREX_ASSERT(AMREX_SPACEDIM==3);
    
    // 
    // First allocate fc
    //
    BoxArray x_ba = cc.boxArray();
    x_ba.surroundingNodes(0);
    fc[0].reset(new MultiFab(x_ba,cc.DistributionMap(),1,nghost, MFInfo(), *ebfactory[lev]));
    fc[0]->setVal(1.0e200);

    BoxArray y_ba = cc.boxArray();
    y_ba.surroundingNodes(1);
    fc[1].reset(new MultiFab(y_ba,cc.DistributionMap(),1,nghost, MFInfo(), *ebfactory[lev]));
    fc[1]->setVal(1.0e200);

    BoxArray z_ba = cc.boxArray();
    z_ba.surroundingNodes(2);
    fc[2].reset(new MultiFab(z_ba,cc.DistributionMap(),1,nghost, MFInfo(), *ebfactory[lev]));
    fc[2]->setVal(1.0e200);

    //
    // Average
    // We do not care about EB because faces in covered regions
    // should never get used so we can set them to whatever values
    // we like
    //
#ifdef _OPENMP
#pragma omp parallel 
#endif
    for (MFIter mfi(*vel[lev],true); mfi.isValid(); ++mfi)
    {
        // Boxes for staggered components
        Box bx = mfi.tilebox();

        average_cc_to_fc( BL_TO_FORTRAN_BOX(bx),
                          BL_TO_FORTRAN_ANYD((*fc[0])[mfi]),
                          BL_TO_FORTRAN_ANYD((*fc[1])[mfi]),
                          BL_TO_FORTRAN_ANYD((*fc[2])[mfi]),
                          BL_TO_FORTRAN_ANYD(cc[mfi]));
    }

    // Set BCs for the fac-centered velocities
	fc[0]->FillBoundary(geom[lev].periodicity());
	fc[1]->FillBoundary(geom[lev].periodicity());
	fc[2]->FillBoundary(geom[lev].periodicity());

	Box domain(geom[lev].Domain());

#ifdef _OPENMP
#pragma omp parallel 
#endif
    for (MFIter mfi(*vel[lev],true); mfi.isValid(); ++mfi)
    {
        // Boxes for staggered components
        Box bx = mfi.tilebox();

     	set_mac_velocity_bcs(bx.loVect(),
     						 bx.hiVect(),
     						 BL_TO_FORTRAN_ANYD((*fc[0])[mfi]),
     						 BL_TO_FORTRAN_ANYD((*fc[1])[mfi]),
     						 BL_TO_FORTRAN_ANYD((*fc[2])[mfi]),
     						 bc_ilo.dataPtr(),
     						 bc_ihi.dataPtr(),
     						 bc_jlo.dataPtr(),
     						 bc_jhi.dataPtr(),
     						 bc_klo.dataPtr(),
     						 bc_khi.dataPtr(),
     						 domain.loVect(),
     						 domain.hiVect(),
     						 &nghost);
    }
    
} 
