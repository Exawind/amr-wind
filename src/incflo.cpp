#include <AMReX_EBAmrUtil.H>
#include <AMReX_EBMultiFabUtil.H>
#include <AMReX_ParmParse.H>

#include <incflo.H>

// Declare and initialise variables
Real stop_time = -1.0;
int max_step = -1;
bool steady_state = false;

int check_int = -1;
int last_chk = -1;
std::string check_file{"chk"};
std::string restart_file{""};

int plot_int = -1;
int last_plt = -1;
std::string plot_file{"plt"};
bool write_eb_surface = false;

int repl_x = 1; int repl_y = 1; int repl_z = 1;
int regrid_int = -1;

void ReadParameters()
{
	{
		ParmParse pp("amr");

		pp.query("stop_time", stop_time);
		pp.query("max_step", max_step);
		pp.query("steady_state", steady_state);

		pp.query("check_file", check_file);
		pp.query("check_int", check_int);
		pp.query("restart", restart_file);

		pp.query("plot_file", plot_file);
		pp.query("plot_int", plot_int);
		pp.query("write_eb_surface", write_eb_surface);

		pp.query("repl_x", repl_x);
		pp.query("repl_y", repl_y);
		pp.query("repl_z", repl_z);
		pp.query("regrid_int", regrid_int);
	}
}

// Initiate vars which cannot be initiated in header
Vector<Real> incflo::gravity(3, 0.);
std::string incflo::load_balance_type = "FixedSize";
std::string incflo::knapsack_weight_type = "RunTimeCosts";

// Define unit vectors for easily convert indeces
amrex::IntVect incflo::e_x(1, 0, 0);
amrex::IntVect incflo::e_y(0, 1, 0);
amrex::IntVect incflo::e_z(0, 0, 1);

int incflo::nlev = 1;

EBSupport incflo::m_eb_support_level = EBSupport::full;


// Constructor
incflo::incflo()
{
    // Read parameters from ParmParse
    ReadParameters();

    // Geometry on all levels has just been defined in the AmrCore constructor

    // No valid BoxArray and DistributionMapping have been defined.
    // But the arrays for them have been resized.

    nlev = maxLevel() + 1;
    istep.resize(nlev, 0);
}

incflo::~incflo(){};

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

            if(!steady_state && regrid_int > -1 && nstep % regrid_int == 0)
                Regrid();

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

void incflo::InitData()
{
	// Initialize internals from ParmParse database
	InitParams();

	// Initialize memory for data-array internals
	ResizeArrays();

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
		IntVect Nrep(repl_x, repl_y, repl_z);
		Restart(restart_file, &nstep, &dt, &time, Nrep);
	}


	// Regrid
	Regrid();

    // Post-initialisation step
	PostInit(dt, time, nstep, restart_flag, stop_time, steady_state);

	// Write out EB sruface
	if(write_eb_surface)
		WriteEBSurface();
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


void incflo::Regrid()
{
	BL_PROFILE_REGION_START("incflo::Regrid()");

    int base_lev = 0;

	if(load_balance_type == "KnapSack")
	{
		AMREX_ALWAYS_ASSERT(fluid_cost[0] != nullptr);

		if(ParallelDescriptor::NProcs() == 1)
			return;
		{
			MultiFab costs(grids[base_lev], dmap[base_lev], 1, 0);
			costs.setVal(0.0);
			costs.plus(*fluid_cost[base_lev], 0, 1, 0);

			DistributionMapping newdm = DistributionMapping::makeKnapSack(costs);

			bool dm_changed = (newdm != dmap[base_lev]);

			SetDistributionMap(base_lev, newdm);

			if(dm_changed)
				RegridArrays(base_lev);

			fluid_cost[base_lev].reset(new MultiFab(grids[base_lev], newdm, 1, 0));
			fluid_cost[base_lev]->setVal(0.0);

			incflo_set_bc0();

			if(ebfactory[base_lev])
			{
				ebfactory[base_lev].reset(new EBFArrayBoxFactory(
					*eb_level_fluid,
					geom[base_lev],
					grids[base_lev],
					dmap[base_lev],
					{m_eb_basic_grow_cells, m_eb_volume_grow_cells, m_eb_full_grow_cells},
					m_eb_support_level));
			}
		}
	}
	BL_PROFILE_REGION_STOP("incflo::Regrid()");
}
