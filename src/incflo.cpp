#include <AMReX_EBAmrUtil.H>
#include <AMReX_EBMultiFabUtil.H>

#include <incflo.H>

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
	// Initialize derived internals
	Init();

	// Either init from scratch or from the checkpoint file
	int restart_flag = 0;
	if(restart_file.empty())
	{
		// NOTE: this also builds ebfactories
		InitLevelData(time);
	}
	else
	{
		restart_flag = 1;

		// NOTE: 1) this also builds ebfactories
        //       2) this can change the grids (during replication)
		Restart(restart_file, &nstep, &dt, &time);
	}

    // Post-initialisation step
	PostInit(dt, time, nstep, restart_flag, stop_time, steady_state);

	// Write out EB sruface
	if(write_eb_surface)
		WriteEBSurface();
}

void incflo::Evolve()
{
    BL_PROFILE("Evolve");
    BL_PROFILE_REGION("Evolve");

	int finish = 0;

	// Initialize prev_dt here; it will be re-defined by call to evolve_fluid
	Real prev_dt = dt;

	// We automatically write checkpoint and plotfiles with the initial data
	//    if plot_int > 0
	if(restart_file.empty() && plot_int > 0)
	{
		incflo_compute_strainrate();
		incflo_compute_vort();
		WritePlotFile(plot_file, nstep, dt, time);
	}

	// We automatically write checkpoint files with the initial data
	//    if check_int > 0
	if(restart_file.empty() && check_int > 0)
	{
		WriteCheckPointFile(check_file, nstep, dt, time);
		last_chk = nstep;
	}

	bool do_not_evolve =
		!steady_state && ((max_step == 0) || ((stop_time >= 0.) && (time > stop_time)) ||
						  ((stop_time <= 0.) && (max_step <= 0)));

    if(!do_not_evolve)
    {
        while(finish == 0)
        {
            Real strt_step = ParallelDescriptor::second();

            Advance(nstep, steady_state, dt, prev_dt, time, stop_time);

            Real end_step = ParallelDescriptor::second() - strt_step;
            ParallelDescriptor::ReduceRealMax(end_step,
                                              ParallelDescriptor::IOProcessorNumber());
            amrex::Print() << "Time per step " << end_step << std::endl;

            if(!steady_state)
            {
                time += prev_dt;
                nstep++;

                if((plot_int > 0) && (nstep % plot_int == 0))
                {
                    incflo_compute_strainrate();
                    incflo_compute_vort();
                    WritePlotFile(plot_file, nstep, dt, time);
                    last_plt = nstep;
                }

                if((check_int > 0) && (nstep % check_int == 0))
                {
                    WriteCheckPointFile(check_file, nstep, dt, time);
                    last_chk = nstep;
                }
            }

            // Mechanism to terminate incflo normally.
            do_not_evolve =
                steady_state || (((stop_time >= 0.) && (time + 0.1 * dt >= stop_time)) ||
                                 (max_step >= 0 && nstep >= max_step));
            if(do_not_evolve)
                finish = 1;
        }
    }

	if(steady_state)
		nstep = 1;

	// Dump plotfile at the final time
	if(check_int > 0 && nstep != last_chk)
		WriteCheckPointFile(check_file, nstep, dt, time);
	if(plot_int > 0 && nstep != last_plt)
    {
		incflo_compute_strainrate();
        incflo_compute_vort();
		WritePlotFile(plot_file, nstep, dt, time);
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
