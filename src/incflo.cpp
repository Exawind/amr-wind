
#include <incflo.H>

// Need this for TagCutCells
#ifdef AMREX_USE_EB
#include <AMReX_EBAmrUtil.H>
#endif

using namespace amrex;

// Constructor
// Note that geometry on all levels has already been defined in the AmrCore constructor,
// which the incflo class inherits from.
incflo::incflo ()
{
    // NOTE: Geometry on all levels has just been defined in the AmrCore
    // constructor. No valid BoxArray and DistributionMapping have been defined.
    // But the arrays for them have been resized.

    // Read inputs file using ParmParse
    ReadParameters();

#ifdef AMREX_USE_EB
    // This is needed before initializing level MultiFabs: ebfactories should
    // not change after the eb-dependent MultiFabs are allocated.
    MakeEBGeometry();
#endif

    // Initialize memory for data-array internals
    ResizeArrays();

    init_bcs();

    set_background_pressure();

    // xxxxx flux registers
}

incflo::~incflo ()
{}

void incflo::InitData ()
{
    BL_PROFILE("incflo::InitData()");

    int restart_flag = 0;
    if(restart_file.empty())
    {
        // This tells the AmrMesh class not to iterate when creating the initial
        // grid hierarchy
        // SetIterateToFalse();

        // This tells the Cluster routine to use the new chopping routine which
        // rejects cuts if they don't improve the efficiency
        SetUseNewChop();

        // This is an AmrCore member function which recursively makes new levels
        // with MakeNewLevelFromScratch.
        InitFromScratch(cur_time);

        if (do_initial_proj) {
            InitialProjection();
        }
        if (initial_iterations > 0) {
            InitialIterations();
        }

        // xxxxx averagedown ???

        // xxxxx if (check_int > 0) { WriteCheckPointFile(); }
    }
    else
    {
        amrex::Abort("xxxxx restart todo");
        // Read starting configuration from chk file.
        ReadCheckpointFile();
    }

    // Post-initialisation step
    // - Initialize diffusive and projection operators
    // - Fill boundaries
    // - Create instance of MAC projection class
    // - Apply initial conditions
    // - Project initial velocity to make divergence free
    // - Perform dummy iterations to find pressure distribution
//    PostInit(restart_flag);

    // Plot initial distribution
    if((plot_int > 0 || plot_per_exact > 0 || plot_per_approx > 0) && !restart_flag)
    {
        amrex::Warning("xxxxx Plotfile todo");
//        WritePlotFile();
        last_plt = 0;
    }
    if(KE_int > 0 && !restart_flag)
    {
        amrex::Abort("xxxxx KE_int todo");
//        amrex::Print() << "Time, Kinetic Energy: " << cur_time << ", " << ComputeKineticEnergy() << std::endl;
    }

#ifdef AMREX_USE_EB
    ParmParse pp("incflo");
    bool write_eb_surface = false;
    pp.query("write_eb_surface", write_eb_surface);
    if (write_eb_surface)
    {
        amrex::Print() << "Writing the geometry to a vtp file.\n" << std::endl;
        amrex::Warning("xxxxx WriteMyEBSurface todo");
  //      WriteMyEBSurface();
    }
#endif
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
                if (nstep % regrid_int == 0)
                {
                    regrid(lev, time);
                    incflo_setup_solvers();
                }
         
            }
         
            if (nstep % regrid_int == 0)
            {
              setup_level_mask();
            }
         
        }*/

        // Advance to time t + dt
        Advance();
        nstep++;
        cur_time += dt;

        amrex::Abort("xxxxx After the first Advance()");

        if (writeNow())
        {
            WritePlotFile();
            last_plt = nstep;
        }

        if(check_int > 0 && (nstep % check_int == 0))
        {
            WriteCheckPointFile();
            last_chk = nstep;
        }
        
        if(KE_int > 0 && (nstep % KE_int == 0))
        {
            amrex::Print() << "Time, Kinetic Energy: " << cur_time << ", " << ComputeKineticEnergy() << std::endl;
        }

        // Mechanism to terminate incflo normally.
        do_not_evolve = (steady_state && SteadyStateReached()) ||
                        ((stop_time > 0. && (cur_time >= stop_time - 1.e-12 * dt)) ||
                         (max_step >= 0 && nstep >= max_step));
    }

    amrex::Abort("xxxxx At the end of Evolve");

	// Output at the final time
    if( check_int > 0                                               && nstep != last_chk) WriteCheckPointFile();
    if( (plot_int > 0 || plot_per_exact > 0 || plot_per_approx > 0) && nstep != last_plt) WritePlotFile();
}

// tag cells for refinement
// overrides the pure virtual function in AmrCore
void incflo::ErrorEst(int lev,
                      TagBoxArray& tags,
                      Real time,
                      int ngrow)
{
    BL_PROFILE("incflo::ErrorEst()");

    const char   tagval = TagBox::SET;
    const char clearval = TagBox::CLEAR;

#ifdef AMREX_USE_EB
    auto const& factory = dynamic_cast<EBFArrayBoxFactory const&>(vel[lev]->Factory());
    auto const& flags = factory.getMultiEBCellFlagFab();
#endif

    const Real* dx      = geom[lev].CellSize();
    const Real* prob_lo = geom[lev].ProbLo();

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(*vel[lev],true); mfi.isValid(); ++mfi)
    {
#ifdef AMREX_USE_EB
        const Box& bx  = mfi.tilebox();
        const auto& flag = flags[mfi];
        const FabType typ = flag.getType(bx);
        if (typ != FabType::covered)
        {
            TagBox&     tagfab  = tags[mfi];

            amrex::Abort("xxxxx TODO: ErrorEst");
#if 0
            // tag cells for refinement
            state_error(BL_TO_FORTRAN_BOX(bx),
                        BL_TO_FORTRAN_ANYD(tagfab),
                        BL_TO_FORTRAN_ANYD((ebfactory[lev]->getVolFrac())[mfi]),
                        &tagval, &clearval,
                        AMREX_ZFILL(dx), AMREX_ZFILL(prob_lo), &time);
#endif
        }
#else
            TagBox&     tagfab  = tags[mfi];

            amrex::Abort("xxxxx TODO: ErrorEst");

            // tag cells for refinement
//          state_error(BL_TO_FORTRAN_BOX(bx),
//                      BL_TO_FORTRAN_ANYD(tagfab),
//                      BL_TO_FORTRAN_ANYD((ebfactory[lev]->getVolFrac())[mfi]),
//                      &tagval, &clearval,
//                      AMREX_ZFILL(dx), AMREX_ZFILL(prob_lo), &time);
#endif
    }

#ifdef AMREX_USE_EB
    refine_cutcells = true;
    // Refine on cut cells
    if (refine_cutcells)
    {
        amrex::TagCutCells(tags, *vel[lev]);
    }
#endif
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

#ifdef AMREX_USE_EB
    m_factory[lev] = makeEBFabFactory(geom[lev], grids[lev], dmap[lev],
                                      {nghost_eb_basic(),
                                       nghost_eb_volume(),
                                       nghost_eb_full()},
                                       EBSupport::full);
#else
    m_factory[lev].reset(new FArrayBoxFactory());
#endif

    m_leveldata[lev].reset(new LevelData(grids[lev], dmap[lev], *m_factory[lev],
                                         ntrac, nghost_state(), nghost_force()));

    t_new[lev] = time;
    t_old[lev] = time - 1.e200;

    prob_init_fluid(lev);
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
    amrex::Abort("xxxxx TODO AverageDownTo");
}

bool
incflo::writeNow()
{
    bool write_now = false;

    if ( plot_int > 0 && (nstep % plot_int == 0) ) 
        write_now = true;

    else if ( plot_per_exact  > 0 && (std::abs(std::remainder(cur_time, plot_per_exact)) < 1.e-12) ) 
        write_now = true;

    else if (plot_per_approx > 0.0)
    {
        // Check to see if we've crossed a plot_per_approx interval by comparing
        // the number of intervals that have elapsed for both the current
        // time and the time at the beginning of this timestep.

        int num_per_old = (cur_time-dt) / plot_per_approx;
        int num_per_new = (cur_time   ) / plot_per_approx;

        // Before using these, however, we must test for the case where we're
        // within machine epsilon of the next interval. In that case, increment
        // the counter, because we have indeed reached the next plot_per_approx interval
        // at this point.

        const Real eps = std::numeric_limits<Real>::epsilon() * 10.0 * std::abs(cur_time);
        const Real next_plot_time = (num_per_old + 1) * plot_per_approx;

        if ((num_per_new == num_per_old) && std::abs(cur_time - next_plot_time) <= eps)
        {
            num_per_new += 1;
        }

        // Similarly, we have to account for the case where the old time is within
        // machine epsilon of the beginning of this interval, so that we don't double
        // count that time threshold -- we already plotted at that time on the last timestep.

        if ((num_per_new != num_per_old) && std::abs((cur_time - dt) - next_plot_time) <= eps)
            num_per_old += 1;

        if (num_per_old != num_per_new)
            write_now = true;
    }

    return write_now;
}

Vector<MultiFab*> incflo::get_velocity_old () noexcept
{
    Vector<MultiFab*> r;
    r.reserve(finest_level+1);
    for (int lev = 0; lev <= finest_level; ++lev) {
        r.push_back(&(m_leveldata[lev]->velocity_o));
    }
    return r;
}

Vector<MultiFab*> incflo::get_velocity_new () noexcept
{
    Vector<MultiFab*> r;
    r.reserve(finest_level+1);
    for (int lev = 0; lev <= finest_level; ++lev) {
        r.push_back(&(m_leveldata[lev]->velocity));
    }
    return r;
}

Vector<MultiFab*> incflo::get_density_old () noexcept
{
    Vector<MultiFab*> r;
    r.reserve(finest_level+1);
    for (int lev = 0; lev <= finest_level; ++lev) {
        r.push_back(&(m_leveldata[lev]->density_o));
    }
    return r;
}

Vector<MultiFab*> incflo::get_density_new () noexcept
{
    Vector<MultiFab*> r;
    r.reserve(finest_level+1);
    for (int lev = 0; lev <= finest_level; ++lev) {
        r.push_back(&(m_leveldata[lev]->density));
    }
    return r;
}

Vector<MultiFab*> incflo::get_tracer_old () noexcept
{
    Vector<MultiFab*> r;
    r.reserve(finest_level+1);
    for (int lev = 0; lev <= finest_level; ++lev) {
        r.push_back(&(m_leveldata[lev]->tracer_o));
    }
    return r;
}

Vector<MultiFab*> incflo::get_tracer_new () noexcept
{
    Vector<MultiFab*> r;
    r.reserve(finest_level+1);
    for (int lev = 0; lev <= finest_level; ++lev) {
        r.push_back(&(m_leveldata[lev]->tracer));
    }
    return r;
}

Vector<MultiFab*> incflo::get_vel_forces () noexcept
{
    Vector<MultiFab*> r;
    r.reserve(finest_level+1);
    for (int lev = 0; lev <= finest_level; ++lev) {
        r.push_back(&(m_leveldata[lev]->vel_forces));
    }
    return r;
}

Vector<MultiFab*> incflo::get_tra_forces () noexcept
{
    Vector<MultiFab*> r;
    r.reserve(finest_level+1);
    for (int lev = 0; lev <= finest_level; ++lev) {
        r.push_back(&(m_leveldata[lev]->tra_forces));
    }
    return r;
}

Vector<MultiFab*> incflo::get_conv_velocity_old () noexcept
{
    Vector<MultiFab*> r;
    r.reserve(finest_level+1);
    for (int lev = 0; lev <= finest_level; ++lev) {
        r.push_back(&(m_leveldata[lev]->conv_velocity_o));
    }
    return r;
}

Vector<MultiFab*> incflo::get_conv_velocity_new () noexcept
{
    Vector<MultiFab*> r;
    r.reserve(finest_level+1);
    for (int lev = 0; lev <= finest_level; ++lev) {
        r.push_back(&(m_leveldata[lev]->conv_velocity));
    }
    return r;
}

Vector<MultiFab*> incflo::get_conv_density_old () noexcept
{
    Vector<MultiFab*> r;
    r.reserve(finest_level+1);
    for (int lev = 0; lev <= finest_level; ++lev) {
        r.push_back(&(m_leveldata[lev]->conv_density_o));
    }
    return r;
}

Vector<MultiFab*> incflo::get_conv_density_new () noexcept
{
    Vector<MultiFab*> r;
    r.reserve(finest_level+1);
    for (int lev = 0; lev <= finest_level; ++lev) {
        r.push_back(&(m_leveldata[lev]->conv_density));
    }
    return r;
}

Vector<MultiFab*> incflo::get_conv_tracer_old () noexcept
{
    Vector<MultiFab*> r;
    r.reserve(finest_level+1);
    for (int lev = 0; lev <= finest_level; ++lev) {
        r.push_back(&(m_leveldata[lev]->conv_tracer_o));
    }
    return r;
}

Vector<MultiFab*> incflo::get_conv_tracer_new () noexcept
{
    Vector<MultiFab*> r;
    r.reserve(finest_level+1);
    for (int lev = 0; lev <= finest_level; ++lev) {
        r.push_back(&(m_leveldata[lev]->conv_tracer));
    }
    return r;
}

Vector<MultiFab const*> incflo::get_velocity_old_const () const noexcept
{
    Vector<MultiFab const*> r;
    r.reserve(finest_level+1);
    for (int lev = 0; lev <= finest_level; ++lev) {
        r.push_back(&(m_leveldata[lev]->velocity_o));
    }
    return r;
}

Vector<MultiFab const*> incflo::get_velocity_new_const () const noexcept
{
    Vector<MultiFab const*> r;
    r.reserve(finest_level+1);
    for (int lev = 0; lev <= finest_level; ++lev) {
        r.push_back(&(m_leveldata[lev]->velocity));
    }
    return r;
}

Vector<MultiFab const*> incflo::get_density_old_const () const noexcept
{
    Vector<MultiFab const*> r;
    r.reserve(finest_level+1);
    for (int lev = 0; lev <= finest_level; ++lev) {
        r.push_back(&(m_leveldata[lev]->density_o));
    }
    return r;
}

Vector<MultiFab const*> incflo::get_density_new_const () const noexcept
{
    Vector<MultiFab const*> r;
    r.reserve(finest_level+1);
    for (int lev = 0; lev <= finest_level; ++lev) {
        r.push_back(&(m_leveldata[lev]->density));
    }
    return r;
}

Vector<MultiFab const*> incflo::get_tracer_old_const () const noexcept
{
    Vector<MultiFab const*> r;
    r.reserve(finest_level+1);
    for (int lev = 0; lev <= finest_level; ++lev) {
        r.push_back(&(m_leveldata[lev]->tracer_o));
    }
    return r;
}

Vector<MultiFab const*> incflo::get_tracer_new_const () const noexcept
{
    Vector<MultiFab const*> r;
    r.reserve(finest_level+1);
    for (int lev = 0; lev <= finest_level; ++lev) {
        r.push_back(&(m_leveldata[lev]->tracer));
    }
    return r;
}

Vector<MultiFab const*> incflo::get_vel_forces_const () const noexcept
{
    Vector<MultiFab const*> r;
    r.reserve(finest_level+1);
    for (int lev = 0; lev <= finest_level; ++lev) {
        r.push_back(&(m_leveldata[lev]->vel_forces));
    }
    return r;
}

Vector<MultiFab const*> incflo::get_tra_forces_const () const noexcept
{
    Vector<MultiFab const*> r;
    r.reserve(finest_level+1);
    for (int lev = 0; lev <= finest_level; ++lev) {
        r.push_back(&(m_leveldata[lev]->tra_forces));
    }
    return r;
}

void incflo::copy_from_new_to_old_velocity (IntVect const& ng)
{
    for (int lev = 0; lev <= finest_level; ++lev) {
        copy_from_new_to_old_velocity(lev, ng);
    }
}

void incflo::copy_from_new_to_old_velocity (int lev, IntVect const& ng)
{
    MultiFab::Copy(m_leveldata[lev]->velocity_o,
                   m_leveldata[lev]->velocity, 0, 0, AMREX_SPACEDIM, ng);
}

void incflo::copy_from_old_to_new_velocity (IntVect const& ng)
{
    for (int lev = 0; lev <= finest_level; ++lev) {
        copy_from_old_to_new_velocity(lev, ng);
    }
}

void incflo::copy_from_old_to_new_velocity (int lev, IntVect const& ng)
{
    MultiFab::Copy(m_leveldata[lev]->velocity,
                   m_leveldata[lev]->velocity_o, 0, 0, AMREX_SPACEDIM, ng);
}

void incflo::copy_from_new_to_old_density (IntVect const& ng)
{
    for (int lev = 0; lev <= finest_level; ++lev) {
        copy_from_new_to_old_density(lev, ng);
    }
}

void incflo::copy_from_new_to_old_density (int lev, IntVect const& ng)
{
    MultiFab::Copy(m_leveldata[lev]->density_o,
                   m_leveldata[lev]->density, 0, 0, 1, ng);
}

void incflo::copy_from_old_to_new_density (IntVect const& ng)
{
    for (int lev = 0; lev <= finest_level; ++lev) {
        copy_from_old_to_new_density(lev, ng);
    }
}

void incflo::copy_from_old_to_new_density (int lev, IntVect const& ng)
{
    MultiFab::Copy(m_leveldata[lev]->density,
                   m_leveldata[lev]->density_o, 0, 0, 1, ng);
}

void incflo::copy_from_new_to_old_tracer (IntVect const& ng)
{
    for (int lev = 0; lev <= finest_level; ++lev) {
        copy_from_new_to_old_tracer(lev, ng);
    }
}

void incflo::copy_from_new_to_old_tracer (int lev, IntVect const& ng)
{
    if (ntrac > 0) {
        MultiFab::Copy(m_leveldata[lev]->tracer_o,
                       m_leveldata[lev]->tracer, 0, 0, ntrac, ng);
    }
}

void incflo::copy_from_old_to_new_tracer (IntVect const& ng)
{
    for (int lev = 0; lev <= finest_level; ++lev) {
        copy_from_old_to_new_tracer(lev, ng);
    }
}

void incflo::copy_from_old_to_new_tracer (int lev, IntVect const& ng)
{
    if (ntrac > 0) {
        MultiFab::Copy(m_leveldata[lev]->tracer,
                       m_leveldata[lev]->tracer_o, 0, 0, ntrac, ng);
    }
}
