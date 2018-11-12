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
    // Initialize the IO variables (pltscalarVars etc)
	InitIOData();

	// Either init from scratch or from the checkpoint file
	int restart_flag = 0;
	if(restart_file.empty())
	{
        // This is an AmrCore member function
        InitFromScratch(t);

        // Set covered coarse cells to be the average of overlying fine cells
        AverageDown();

		// NOTE: this also builds ebfactories
		InitLevelData();
	}
	else
	{
		// NOTE: 1) this also builds ebfactories
		ReadCheckpointFile();
		restart_flag = 1;
	}
    
    // Post-initialisation step 
    // - Set BC types 
    // - Fill boundaries
    // - Create instance of MAC projection class
    // - Apply initial conditions
	PostInit(restart_flag);
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

// Set covered coarse cells to be the average of overlying fine cells
// TODO: EB_average_down() does not seem to have support for nodal data -- check with pressure
// TODO: Move somewhere else, for example setup/incflo_arrays.cpp
// TODO: Possible to loop over variables like: "for mf in list of mfs"? 
void incflo::AverageDown()
{
    for (int lev = finest_level - 1; lev >= 0; --lev)
    {
        IntVect rr = refRatio(lev);
        amrex::EB_average_down(*ro[lev+1],          *ro[lev],           0, 1, rr);
        amrex::EB_average_down(*p0[lev+1],          *p0[lev],           0, 1, rr);
        amrex::EB_average_down(*p[lev+1],           *p[lev],            0, 1, rr);
        amrex::EB_average_down(*gp0[lev+1],         *gp0[lev],          0, 3, rr);
        amrex::EB_average_down(*gp[lev+1],          *gp[lev],           0, 3, rr);
        amrex::EB_average_down(*eta[lev+1],         *eta[lev],          0, 1, rr);
        amrex::EB_average_down(*vel[lev+1],         *vel[lev],          0, 3, rr);
        amrex::EB_average_down(*vel_o[lev+1],       *vel_o[lev],        0, 3, rr);
        amrex::EB_average_down(*strainrate[lev+1],  *strainrate[lev],   0, 1, rr);
        amrex::EB_average_down(*vort[lev+1],        *vort[lev],         0, 1, rr);
        amrex::EB_average_down(*divu[lev+1],        *divu[lev],         0, 1, rr);
        amrex::EB_average_down(*phi[lev+1],         *phi[lev],          0, 1, rr);
        amrex::EB_average_down(*phi_diff[lev+1],    *phi_diff[lev],     0, 1, rr);
        amrex::EB_average_down(*rhs_diff[lev+1],    *rhs_diff[lev],     0, 1, rr);
        amrex::EB_average_down(*xslopes[lev+1],     *xslopes[lev],      0, 3, rr);
        amrex::EB_average_down(*yslopes[lev+1],     *yslopes[lev],      0, 3, rr);
        amrex::EB_average_down(*zslopes[lev+1],     *zslopes[lev],      0, 3, rr);
        for (int i = 0; i < 3; i++)
        {
            amrex::EB_average_down(*bcoeff[lev+1][i],       *bcoeff[lev][i],        0, 1, rr);
            amrex::EB_average_down(*bcoeff_diff[lev+1][i],  *bcoeff_diff[lev][i],   0, 1, rr);
        }
        amrex::EB_average_down(*m_u_mac[lev+1], *m_u_mac[lev], 0, 3, rr);
        amrex::EB_average_down(*m_v_mac[lev+1], *m_v_mac[lev], 0, 3, rr);
        amrex::EB_average_down(*m_w_mac[lev+1], *m_w_mac[lev], 0, 3, rr);
    }
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
