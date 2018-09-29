#include <AMReX_ParmParse.H>
#include <AMReX_BC_TYPES.H>
#include <AMReX_Box.H>

#include <incflo_level.H>
#include <eb_F.H>
#include <setup_F.H>

void incflo_level::InitParams()
{
	{
		ParmParse pp("incflo");

        // Verbosity
        pp.query("verbose", verbose);

        // Initial density
		pp.queryarr("gravity", gravity, 0, 3);

        pp.query("ro_0", ro_0);
        AMREX_ASSERT(ro_0 >= 0.0);

        pp.query("mu_0", mu_0);
        AMREX_ASSERT(mu_0 >= 0.0);

		// Options to control time stepping
		pp.query("cfl", cfl);
		pp.query("fixed_dt", fixed_dt);

		// Tolerance to check for steady state
		pp.query("steady_state_tol", steady_state_tol);

		// Option to control MGML behavior
		pp.query("mg_verbose", mg_verbose);
		pp.query("mg_cg_verbose", mg_cg_verbose);
		pp.query("mg_max_iter", mg_max_iter);
		pp.query("mg_cg_maxiter", mg_cg_maxiter);
		pp.query("mg_max_fmg_iter", mg_max_fmg_iter);
		pp.query("mg_rtol", mg_rtol);
		pp.query("mg_atol", mg_atol);

        // Default bottom solver is bicgstab, but alternatives are "smoother" or "hypre"
        bottom_solver_type = "bicgstab";
        pp.query( "bottom_solver_type",  bottom_solver_type );

		// Should we use explicit vs implicit diffusion
		pp.query("explicit_diffusion", explicit_diffusion);

		// If we are using the MAC-projected velocities, how should we discretize the
		//  ugradu term.
		pp.query("ugradu_type", ugradu_type);

		// The default type is "FixedSize" but we can over-write that in the inputs file
		//  with "KDTree" or "KnapSack"
		pp.query("load_balance_type", load_balance_type);
		pp.query("knapsack_weight_type", knapsack_weight_type);

		AMREX_ASSERT(load_balance_type == "FixedSize" || load_balance_type == "KnapSack");
		AMREX_ASSERT(knapsack_weight_type == "RunTimeCosts");
        
        // Loads constants given at runtime `inputs` file into the Fortran module "constant"
        incflo_get_data(gravity.dataPtr(), &ro_0, &mu_0); 
	}
}

void incflo_level::Init(int lev, Real time)
{
	BL_ASSERT(max_level == 0);

	InitIOData();

	// Define coarse level BoxArray and DistributionMap
	finest_level = 0;

	const BoxArray& ba = MakeBaseGrids();
	DistributionMapping dm(ba, ParallelDescriptor::NProcs());

	MakeNewLevelFromScratch(0, time, ba, dm);

	int cyc_x = 0, cyc_y = 0, cyc_z = 0;
	if(geom[lev].isPeriodic(0))
		cyc_x = 1;
	if(geom[lev].isPeriodic(1))
		cyc_y = 1;
	if(geom[lev].isPeriodic(2))
		cyc_z = 1;

	incflo_set_cyclic(&cyc_x, &cyc_y, &cyc_z);

	incflo_set_bc_type(lev);

	// Allocate container for eb-normals
	dummy = std::unique_ptr<MultiFab>(new MultiFab);

	// Create MAC projection object
	mac_projection.reset(new MacProjection(this, nghost, &ebfactory));
	mac_projection->set_bcs(&bc_ilo, &bc_ihi, &bc_jlo, &bc_jhi, &bc_klo, &bc_khi);
}

BoxArray incflo_level::MakeBaseGrids() const
{
	BoxArray ba(geom[0].Domain());

	ba.maxSize(max_grid_size[0]);

	// We only call ChopGrids if dividing up the grid using max_grid_size didn't
	//    create enough grids to have at least one grid per processor.
	// This option is controlled by "refine_grid_layout" which defaults to true.

	if(refine_grid_layout && ba.size() < ParallelDescriptor::NProcs() &&
	   (load_balance_type == "FixedSize" || load_balance_type == "KnapSack"))
	{
		ChopGrids(geom[0].Domain(), ba, ParallelDescriptor::NProcs());
	}

	if(ba == grids[0])
	{
		ba = grids[0]; // to avoid dupliates
	}
	amrex::Print() << "In MakeBaseGrids: BA HAS " << ba.size() << " GRIDS " << std::endl;
	return ba;
}

void incflo_level::ChopGrids(const Box& domain, BoxArray& ba, int target_size) const
{
	if(ParallelDescriptor::IOProcessor())
		amrex::Warning(
			"Using max_grid_size only did not make enough grids for the number of processors");

	// Here we hard-wire the maximum number of times we divide the boxes.
	int n = 10;

	// Here we hard-wire the minimum size in any one direction the boxes can be
	int min_grid_size = 4;

	IntVect chunk(domain.length(0), domain.length(1), domain.length(2));

	int j;
	for(int cnt = 1; cnt <= n; ++cnt)
	{
		if(chunk[0] >= chunk[1] && chunk[0] >= chunk[2])
		{
			j = 0;
		}
		else if(chunk[1] >= chunk[0] && chunk[1] >= chunk[2])
		{
			j = 1;
		}
		else if(chunk[2] >= chunk[0] && chunk[2] >= chunk[1])
		{
			j = 2;
		}
		chunk[j] /= 2;

		if(chunk[j] >= min_grid_size)
		{
			ba.maxSize(chunk);
		}
		else
		{
			// chunk[j] was the biggest chunk -- if this is too small then we're done
			if(ParallelDescriptor::IOProcessor())
				amrex::Warning(
					"ChopGrids was unable to make enough grids for the number of processors");
			return;
		}

		// Test if we now have enough grids
		if(ba.size() >= target_size)
			return;
	}
}

void incflo_level::MakeNewLevelFromScratch(int lev,
										   Real time,
										   const BoxArray& new_grids,
										   const DistributionMapping& new_dmap)
{
	SetBoxArray(lev, new_grids);
	SetDistributionMap(lev, new_dmap);

	MakeBCArrays();
}

void incflo_level::ReMakeNewLevelFromScratch(int lev,
											 const BoxArray& new_grids,
											 const DistributionMapping& new_dmap)
{
	SetBoxArray(lev, new_grids);
	SetDistributionMap(lev, new_dmap);

	MakeBCArrays();

	// We need to re-fill these arrays for the larger domain (after replication).
	incflo_set_bc_type(lev);
}

void incflo_level::InitLevelData(int lev, Real time)
{
	// This needs is needed before initializing level MultiFabs: ebfactories should
	// not change after the eb-dependent MultiFabs are allocated.
	make_eb_geometry(lev);

	// Allocate the fluid data, NOTE: this depends on the ebfactories.
	AllocateArrays(lev);
}

void incflo_level::PostInit(
	int lev, Real& dt, Real time, int nstep, int restart_flag, Real stop_time, int steady_state)
{
    // Initial fluid arrays: pressure, velocity, density, viscosity
    incflo_init_fluid(lev, restart_flag, time, dt, stop_time, steady_state);
}

void incflo_level::MakeBCArrays()
{
	// Define and allocate the integer MultiFab that is the outside adjacent cells of the problem domain.
	Box domainx(geom[0].Domain());
	domainx.grow(1, nghost);
	domainx.grow(2, nghost);
	Box box_ilo = amrex::adjCellLo(domainx, 0, 1);
	Box box_ihi = amrex::adjCellHi(domainx, 0, 1);

	Box domainy(geom[0].Domain());
	domainy.grow(0, nghost);
	domainy.grow(2, nghost);
	Box box_jlo = amrex::adjCellLo(domainy, 1, 1);
	Box box_jhi = amrex::adjCellHi(domainy, 1, 1);

	Box domainz(geom[0].Domain());
	domainz.grow(0, nghost);
	domainz.grow(1, nghost);
	Box box_klo = amrex::adjCellLo(domainz, 2, 1);
	Box box_khi = amrex::adjCellHi(domainz, 2, 1);

	// Note that each of these is a single IArrayBox so every process has a copy of them
	bc_ilo.resize(box_ilo, 2);
	bc_ihi.resize(box_ihi, 2);
	bc_jlo.resize(box_jlo, 2);
	bc_jhi.resize(box_jhi, 2);
	bc_klo.resize(box_klo, 2);
	bc_khi.resize(box_khi, 2);
}

void incflo_level::incflo_init_fluid(int lev, int is_restarting, Real time, Real& dt, 
                                     Real stop_time, int steady_state)
{
	Box domain(geom[lev].Domain());

	Real dx = geom[lev].CellSize(0);
	Real dy = geom[lev].CellSize(1);
	Real dz = geom[lev].CellSize(2);

	Real xlen = geom[lev].ProbHi(0) - geom[lev].ProbLo(0);
	Real ylen = geom[lev].ProbHi(1) - geom[lev].ProbLo(1);
	Real zlen = geom[lev].ProbHi(2) - geom[lev].ProbLo(2);

	// Here we set bc values for p and u,v,w before the IC's are set
	incflo_set_bc0(lev);

	// We deliberately don't tile this loop since we will be looping
	//    over bc's on faces and it makes more sense to do this one grid at a time
	for(MFIter mfi(*ro[lev], false); mfi.isValid(); ++mfi)
	{

		const Box& bx = mfi.validbox();
		const Box& sbx = (*ro[lev])[mfi].box();

		if(is_restarting)
		{

			init_fluid_restart(sbx.loVect(),
							   sbx.hiVect(),
							   bx.loVect(),
							   bx.hiVect(),
							   (*mu[lev])[mfi].dataPtr(),
							   (*lambda[lev])[mfi].dataPtr());
		}
		else
		{

			init_fluid(sbx.loVect(),
                       sbx.hiVect(),
					   bx.loVect(),
					   bx.hiVect(),
					   domain.loVect(),
					   domain.hiVect(),
					   (*ro[lev])[mfi].dataPtr(),
					   (*p[lev])[mfi].dataPtr(),
					   (*vel[lev])[mfi].dataPtr(),
					   (*mu[lev])[mfi].dataPtr(),
					   (*lambda[lev])[mfi].dataPtr(),
                       &dx, &dy, &dz,
                       &xlen, &ylen, &zlen);
		}
	}

	incflo_set_p0(lev);

	// Here we re-set the bc values for p and u,v,w just in case init_fluid
	//      over-wrote some of the bc values with ic values
	incflo_set_bc0(lev);

	// We deliberately don't tile this loop since we will be looping
	//    over bc's on faces and it makes more sense to do this one grid at a time
	if(!is_restarting)
	{

		for(MFIter mfi(*ro[lev]); mfi.isValid(); ++mfi)
		{

			const Box& sbx = (*ro[lev])[mfi].box();

			zero_wall_norm_vel(sbx.loVect(),
							   sbx.hiVect(),
							   (*vel[lev])[mfi].dataPtr(),
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

	if(!nodal_pressure)
		incflo_extrap_pressure(lev, p0[lev]);

	fill_mf_bc(lev, *ro[lev]);

	vel[lev]->FillBoundary(geom[lev].periodicity());

	fill_mf_bc(lev, *mu[lev]);
	fill_mf_bc(lev, *lambda[lev]);

	if(is_restarting == 1)
	{
		if(!nodal_pressure)
			incflo_extrap_pressure(lev, p[lev]);
	}
	else
	{
		// Here initialize dt to -1 so that we don't check new dt against a previous value
		dt = -1.;
		incflo_compute_dt(lev, time, stop_time, steady_state, dt);

		incflo_set_scalar_bcs(lev);
		incflo_project_velocity(lev);
		incflo_initial_iterations(lev, dt, stop_time, steady_state);
	}
}

void incflo_level::incflo_set_bc0(int lev)
{
	Box domain(geom[lev].Domain());

	// Don't tile this -- at least for now
	for(MFIter mfi(*ro[lev]); mfi.isValid(); ++mfi)
	{
		const Box& sbx = (*ro[lev])[mfi].box();

		set_bc0(sbx.loVect(),
				sbx.hiVect(),
				(*ro[lev])[mfi].dataPtr(),
				(*mu[lev])[mfi].dataPtr(),
				(*lambda[lev])[mfi].dataPtr(),
				bc_ilo.dataPtr(),
				bc_ihi.dataPtr(),
				bc_jlo.dataPtr(),
				bc_jhi.dataPtr(),
				bc_klo.dataPtr(),
				bc_khi.dataPtr(),
				domain.loVect(),
				domain.hiVect(),
				&nghost,
				&nodal_pressure);
	}

	if(!nodal_pressure)
		fill_mf_bc(lev, *p[lev]);

	fill_mf_bc(lev, *ro[lev]);

	// Put velocity Dirichlet bc's on faces
	int extrap_dir_bcs = 0;
	incflo_set_velocity_bcs(lev, extrap_dir_bcs);

	vel[lev]->FillBoundary(geom[lev].periodicity());
}

void incflo_level::incflo_set_p0(int lev)
{
	Real xlen = geom[lev].ProbHi(0) - geom[lev].ProbLo(0);
	Real ylen = geom[lev].ProbHi(1) - geom[lev].ProbLo(1);
	Real zlen = geom[lev].ProbHi(2) - geom[lev].ProbLo(2);

	Real dx = geom[lev].CellSize(0);
	Real dy = geom[lev].CellSize(1);
	Real dz = geom[lev].CellSize(2);

	Box domain(geom[lev].Domain());

	int delp_dir;
	set_delp_dir(&delp_dir);

	// We deliberately don't tile this loop since we will be looping
	//    over bc's on faces and it makes more sense to do this one grid at a time
	for(MFIter mfi(*ro[lev], false); mfi.isValid(); ++mfi)
	{

		const Box& bx = mfi.validbox();

        set_p0(bx.loVect(), bx.hiVect(),
               domain.loVect(), domain.hiVect(),
			   BL_TO_FORTRAN_ANYD((*p0[lev])[mfi]),
			   BL_TO_FORTRAN_ANYD((*gp0[lev])[mfi]),
               &dx, &dy, &dz,
               &xlen, &ylen, &zlen,
			   &delp_dir,
               bc_ilo.dataPtr(), bc_ihi.dataPtr(),
               bc_jlo.dataPtr(), bc_jhi.dataPtr(),
               bc_klo.dataPtr(), bc_khi.dataPtr(),
			   &nghost,
			   &nodal_pressure);
	}

	// Here we set a separate periodicity flag for p0 because when we use
	// pressure drop (delp) boundary conditions we fill all variables *except* p0
	// periodically

	IntVect press_per =
		IntVect(geom[lev].isPeriodic(0), geom[lev].isPeriodic(1), geom[lev].isPeriodic(2));

	if(delp_dir > -1)
		press_per[delp_dir] = 0;
	p0_periodicity = Periodicity(press_per);

	p0[lev]->FillBoundary(p0_periodicity);
	 gp0[lev]->FillBoundary(p0_periodicity);
}

// 
// Perform initial pressure iterations 
//
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

