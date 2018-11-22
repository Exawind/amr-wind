#include <AMReX_EBAmrUtil.H>
#include <AMReX_EBMultiFabUtil.H>

#include <incflo.H>
#include <derive_F.H>

// Constructor
// Note that geometry on all levels has already been defined in the AmrCore constructor,
// which the incflo class inherits from.
incflo::incflo()
{
    // Read inputs file using ParmParse
    ReadParameters();

	// Initialize memory for data-array internals
	ResizeArrays();

    // This needs is needed before initializing level MultiFabs: ebfactories should
    // not change after the eb-dependent MultiFabs are allocated.
    MakeEBGeometry();
    if(incflo_verbose > 0) WriteEBSurface();
}

incflo::~incflo(){};

void incflo::InitData()
{
    BL_PROFILE("incflo::InitData()");

	// Either init from scratch or from the checkpoint file
	int restart_flag = 0;
	if(restart_file.empty())
	{
        // This is an AmrCore member function. 
        // Importantly, it calls MakeNewLevelFromScratch():
        // - Set BA and DM 
        // - Allocate arrays for level
        InitFromScratch(cur_time);
	}
	else
	{
        // Read starting configuration from chk file. 
        // Importantly, it calls MakeNewLevelFromScratch():
        // - Set BA and DM 
        // - Allocate arrays for level
		ReadCheckpointFile();
		restart_flag = 1;
	}

    // Post-initialisation step
    // - Set BC types
    // - Fill boundaries
    // - Create instance of MAC projection class
    // - Apply initial conditions
    // - Project initial velocity to make divergence free
    // - Perform dummy iterations to find pressure distribution
	PostInit(restart_flag);

    // Plot initial distribution
    if(plot_int > 0 and !restart_flag) 
    {
        UpdateDerivedQuantities();
        WritePlotFile();
        last_plt = 0;
    }
}

void incflo::Evolve()
{
    BL_PROFILE("incflo::Evolve()");

    bool do_not_evolve = ((max_step == 0) || ((stop_time >= 0.) && (cur_time > stop_time)) ||
						  ((stop_time <= 0.) && (max_step <= 0))) && !steady_state;

    while(!do_not_evolve)
    {
        // TODO: Necessary for dynamic meshing
        /* if (regrid_int > 0) 
        {
            // Make sure we don't regrid on max_level
            for (int lev = 0; lev < max_level; ++lev)
            {
                // regrid is a member function of AmrCore
                if (nstep % regrid_int == 0) regrid(lev, time);
            }
        }*/

        // Advance to time t + dt
        Advance();
        nstep++;
        cur_time += dt;

        // Write plot and checkpoint files
        if((plot_int > 0) && (nstep % plot_int == 0))
        {
            UpdateDerivedQuantities();
            WritePlotFile();
            last_plt = nstep;
        }
        if((check_int > 0) && (nstep % check_int == 0))
        {
            WriteCheckPointFile();
            last_chk = nstep;
        }

        // Mechanism to terminate incflo normally.
        do_not_evolve = (steady_state && SteadyStateReached()) ||
                        (((stop_time > 0.) && (cur_time >= stop_time - 1.e-6 * dt)) ||
                         (max_step >= 0 && nstep >= max_step));
    }

	// Output at the final time
    if(check_int > 0 && nstep != last_chk) WriteCheckPointFile();
    if(plot_int > 0 && nstep != last_plt)
    {
        UpdateDerivedQuantities();
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

    const int clearval = TagBox::CLEAR;
    const int   tagval = TagBox::SET;

    const Real* dx      = geom[lev].CellSize();
    const Real* prob_lo = geom[lev].ProbLo();

    Vector<int>  itags;

    for (MFIter mfi(*ro[lev],true); mfi.isValid(); ++mfi)
    {
        const Box& tilebox  = mfi.tilebox();

        TagBox&     tagfab  = tags[mfi];
        
        // We cannot pass tagfab to Fortran becuase it is BaseFab<char>.
        // So we are going to get a temporary integer array.
            // set itags initially to 'untagged' everywhere
            // we define itags over the tilebox region
        tagfab.get_itags(itags, tilebox);
        
            // data pointer and index space
        int*        tptr    = itags.dataPtr();
        const int*  tlo     = tilebox.loVect();
        const int*  thi     = tilebox.hiVect();

            // tag cells for refinement
        state_error(tptr,  AMREX_ARLIM_3D(tlo), AMREX_ARLIM_3D(thi),
            BL_TO_FORTRAN_3D((*ro[lev])[mfi]),
            &tagval, &clearval, 
            AMREX_ARLIM_3D(tilebox.loVect()), AMREX_ARLIM_3D(tilebox.hiVect()), 
            AMREX_ZFILL(dx), AMREX_ZFILL(prob_lo), &time);
        //
        // Now update the tags in the TagBox in the tilebox region
            // to be equal to itags
        //
        tagfab.tags_and_untags(itags, tilebox);
    }
    
    /* TODO: This is what we want to refine on, but it gives segfault like this
     * // Refine on cut cells
    if (refine_cutcells) 
    {
        amrex::TagCutCells(tags, *ro[lev]);
    }*/
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

    if(incflo_verbose > 0)
    { 
        amrex::Print() << "Making new level " << lev << std::endl; 
        amrex::Print() << "with BoxArray " << new_grids << std::endl; 
    }

	SetBoxArray(lev, new_grids);
	SetDistributionMap(lev, new_dmap);

    if(lev == 0) MakeBCArrays();

	// Allocate the fluid data, NOTE: this depends on the ebfactories.
    AllocateArrays(lev);
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
// TODO: Move somewhere else, for example setup/incflo_arrays.cpp
void incflo::AverageDown()
{
    BL_PROFILE("incflo::AverageDown()");

    for (int lev = finest_level - 1; lev >= 0; --lev)
    {
        AverageDownTo(lev);
    }
}

void incflo::AverageDownTo(int crse_lev)
{
    BL_PROFILE("incflo::AverageDownTo()");

    IntVect rr = refRatio(crse_lev);
    amrex::EB_average_down(*ro[crse_lev+1],         *ro[crse_lev],         0, 1, rr);
    amrex::EB_average_down(*p0[crse_lev+1],         *p0[crse_lev],         0, 1, rr);
    amrex::EB_average_down(*p[crse_lev+1],          *p[crse_lev],          0, 1, rr);
    amrex::EB_average_down(*eta[crse_lev+1],        *eta[crse_lev],        0, 1, rr);
    amrex::EB_average_down(*strainrate[crse_lev+1], *strainrate[crse_lev], 0, 1, rr);
    amrex::EB_average_down(*vort[crse_lev+1],       *vort[crse_lev],       0, 1, rr);
    amrex::EB_average_down(*gp0[crse_lev+1],        *gp0[crse_lev],        0, 3, rr);
    amrex::EB_average_down(*gp[crse_lev+1],         *gp[crse_lev],         0, 3, rr);
    amrex::EB_average_down(*vel[crse_lev+1],        *vel[crse_lev],        0, 3, rr);
}
