#include <AMReX_EBAmrUtil.H>
#include <AMReX_EBMultiFabUtil.H>

#include <incflo.H>
#include <boundary_conditions_F.H>

// Constructor
// Note that geometry on all levels has already been defined in the AmrCore constructor, 
// which the incflo class inherits from. 
incflo::incflo()
{
    // Read inputs file using ParmParse
    ReadParameters();
    
	// Initialize memory for data-array internals
	ResizeArrays();
}

incflo::~incflo(){};

void incflo::InitData()
{
    // Initialize the IO variables (pltscalarVars etc)
	InitIOData();

	// Either init from scratch or from the checkpoint file
	int restart_flag = 0;
	if(restart_file.empty())
	{
        // These are AmrCore member functions
        InitFromScratch(t);

        // TODO: Implement
        // AverageDown();

		// NOTE: this also builds ebfactories
		InitLevelData();
	}
	else
	{
		// NOTE: 1) this also builds ebfactories
        //       2) this can change the grids (during replication)
		ReadCheckpointFile();
		restart_flag = 1;
	}
    
    // Set BC-types (cyclic only at level 0)
    int cyc_x = 0, cyc_y = 0, cyc_z = 0;
    if(geom[0].isPeriodic(0)) cyc_x = 1;
    if(geom[0].isPeriodic(1)) cyc_y = 1;
    if(geom[0].isPeriodic(2)) cyc_z = 1;
    incflo_set_cyclic(&cyc_x, &cyc_y, &cyc_z);

    for(int lev = 0; lev <= max_level; lev++)
    {
        incflo_set_bc_type(lev);
    }

    // Fill boundaries
	for(int lev = 0; lev <= finest_level; ++lev)
	{
        if(!nodal_pressure) fill_mf_bc(lev, *p[lev]);
		fill_mf_bc(lev, *ro[lev]);
		fill_mf_bc(lev, *eta[lev]);

		// Fill the bc's just in case
		vel[lev]->FillBoundary(geom[lev].periodicity());
		vel_o[lev]->FillBoundary(geom[lev].periodicity());
	}

    // TODO: Put this into PostInit or somethign
    // Create MAC projection object
    mac_projection.reset(new MacProjection(this, nghost, &ebfactory));
    mac_projection->set_bcs(bc_ilo, bc_ihi, bc_jlo, bc_jhi, bc_klo, bc_khi);

    // Post-initialisation step (TODO: this is currrently just init_fluid)
	PostInit(restart_flag);

	// Write out EB sruface
    if(write_eb_surface) WriteEBSurface();
}

void incflo::Evolve()
{
    BL_PROFILE("Evolve");
    BL_PROFILE_REGION("Evolve");

	int finish = 0;

	// We automatically write checkpoint and plotfiles with the initial data
	//    if plot_int > 0
	if(restart_file.empty() && plot_int > 0)
	{
		incflo_compute_strainrate();
		incflo_compute_vort();
		WritePlotFile();
	}

	// We automatically write checkpoint files with the initial data
	//    if check_int > 0
	if(restart_file.empty() && check_int > 0)
	{
		WriteCheckPointFile();
		last_chk = nstep;
	}

	bool do_not_evolve =
		!steady_state && ((max_step == 0) || ((stop_time >= 0.) && (t > stop_time)) ||
						  ((stop_time <= 0.) && (max_step <= 0)));

    if(!do_not_evolve)
    {
        while(finish == 0)
        {
            Real strt_step = ParallelDescriptor::second();

            Advance();

            Real end_step = ParallelDescriptor::second() - strt_step;
            ParallelDescriptor::ReduceRealMax(end_step,
                                              ParallelDescriptor::IOProcessorNumber());
            amrex::Print() << "Time per step " << end_step << std::endl;

            if(!steady_state)
            {
                t += dt;
                nstep++;

                if((plot_int > 0) && (nstep % plot_int == 0))
                {
                    incflo_compute_strainrate();
                    incflo_compute_vort();
                    WritePlotFile();
                    last_plt = nstep;
                }

                if((check_int > 0) && (nstep % check_int == 0))
                {
                    WriteCheckPointFile();
                    last_chk = nstep;
                }
            }

            // Mechanism to terminate incflo normally.
            do_not_evolve =
                steady_state || (((stop_time >= 0.) && fabs(t - stop_time) < 0.01 * dt) ||
                                 (max_step >= 0 && nstep >= max_step));
            if(do_not_evolve)
                finish = 1;
        }
    }

	if(steady_state)
		nstep = 1;

	// Dump plotfile at the final time
	if(check_int > 0 && nstep != last_chk)
		WriteCheckPointFile();
	if(plot_int > 0 && nstep != last_plt)
    {
		incflo_compute_strainrate();
        incflo_compute_vort();
		WritePlotFile();
    }
}

// tag all cells for refinement
// overrides the pure virtual function in AmrCore
void incflo::ErrorEst(int lev,
                      TagBoxArray& tags,
                      Real time,
                      int ngrow)
{
    BL_PROFILE("incflo::ErrorEst()");

    // Refine on cut cells
    if (refine_cutcells) amrex::TagCutCells(tags, *ro[lev]);
}

// Make a new level from scratch using provided BoxArray and DistributionMapping.
// Only used during initialization.
// overrides the pure virtual function in AmrCore
void incflo::MakeNewLevelFromScratch(int lev,
                                     Real time,
                                     const BoxArray& new_grids,
                                     const DistributionMapping& new_dmap)
{
    BL_PROFILE("incflo::MakeNewLevelFromScratch()");

	SetBoxArray(lev, new_grids);
	SetDistributionMap(lev, new_dmap);

    if(lev == 0) MakeBCArrays();
}

// Make a new level using provided BoxArray and DistributionMapping and
// fill with interpolated coarse level data.
// overrides the pure virtual function in AmrCore
void incflo::MakeNewLevelFromCoarse (int lev,
                                     Real time,
                                     const BoxArray& ba,
                                     const DistributionMapping& dm)
{
    BL_PROFILE("incflo::MakeNewLevelFromCoarse()");

    amrex::Print() << "ABORT: incflo::MakeNewLevelFromCoarse() not yet implemented. " << std::endl;
    amrex::Abort();
}

// Remake an existing level using provided BoxArray and DistributionMapping and
// fill with existing fine and coarse data.
// overrides the pure virtual function in AmrCore
void incflo::RemakeLevel (int lev, Real time, const BoxArray& ba,
			 const DistributionMapping& dm)
{
    BL_PROFILE("incflo::RemakeLevel()");

    amrex::Print() << "ABORT: incflo::RemakeLevel() not yet implemented. " << std::endl;
    amrex::Abort();
}

// Delete level data
// overrides the pure virtual function in AmrCore
void incflo::ClearLevel (int lev)
{
    BL_PROFILE("incflo::ClearLevel()");

    amrex::Print() << "ABORT: incflo::ClearLevel() not yet implemented. " << std::endl;
    amrex::Abort();
}

//
// Subroutine to compute norm0 of EB multifab
//
Real incflo::incflo_norm0(const Vector<std::unique_ptr<MultiFab>>& mf, int lev, int comp)
{
	MultiFab mf_tmp(mf[lev]->boxArray(),
					mf[lev]->DistributionMap(),
                    mf[lev]->nComp(),
                    0, MFInfo(), *ebfactory[lev]);

	MultiFab::Copy(mf_tmp, *mf[lev], comp, comp, 1, 0);
	EB_set_covered(mf_tmp, 0.0);

	return mf_tmp.norm0(comp);
}

Real incflo::incflo_norm0(MultiFab& mf, int lev, int comp)
{
    MultiFab mf_tmp(mf.boxArray(),
                    mf.DistributionMap(),
                    mf.nComp(),
                    0, MFInfo(), *ebfactory[lev]);

	MultiFab::Copy(mf_tmp, mf, comp, comp, 1, 0);
	EB_set_covered(mf_tmp, 0.0);

	return mf_tmp.norm0(comp);
}

//
// Subroutine to compute norm1 of EB multifab
//
Real incflo::incflo_norm1(const Vector<std::unique_ptr<MultiFab>>& mf, int lev, int comp)
{
	MultiFab mf_tmp(mf[lev]->boxArray(),
					mf[lev]->DistributionMap(),
					mf[lev]->nComp(),
                    0, MFInfo(), *ebfactory[lev]);

	MultiFab::Copy(mf_tmp, *mf[lev], comp, comp, 1, 0);
	EB_set_covered(mf_tmp, 0.0);

	return mf_tmp.norm1(comp, geom[lev].periodicity());
}

Real incflo::incflo_norm1(MultiFab& mf, int lev, int comp)
{
	MultiFab mf_tmp(mf.boxArray(),
                    mf.DistributionMap(),
                    mf.nComp(),
                    0, MFInfo(), *ebfactory[lev]);

	MultiFab::Copy(mf_tmp, mf, comp, comp, 1, 0);
	EB_set_covered(mf_tmp, 0.0);

	return mf_tmp.norm1(comp, geom[lev].periodicity());
}

//
// Print the maximum values of the velocity components
//
void incflo::incflo_print_max_vel(int lev)
{
	amrex::Print() << "max(abs(u/v/w/p))  = "
                   << incflo_norm0(vel, lev, 0) << "  "
				   << incflo_norm0(vel, lev, 1) << "  "
                   << incflo_norm0(vel, lev, 2) << "  "
				   << incflo_norm0(p, lev, 0) << "  " << std::endl;
}

void incflo::check_for_nans(int lev)
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
