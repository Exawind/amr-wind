
#include <incflo.H>

// Need this for TagCutCells
#ifdef AMREX_USE_EB
#include <AMReX_EBAmrUtil.H>
#endif

using namespace amrex;

incflo::incflo ()
{
    // NOTE: Geometry on all levels has just been defined in the AmrCore
    // constructor. No valid BoxArray and DistributionMapping have been defined.
    // But the arrays for them have been resized.

    // Read inputs file using ParmParse
    ReadParameters();

#ifdef AMREX_USE_EB
    // This is needed before initializing level MultiFab
    MakeEBGeometry();
#endif

    // Initialize memory for data-array internals
    ResizeArrays();

    init_bcs();

    init_advection();

    set_background_pressure();
}

incflo::~incflo ()
{}

void incflo::InitData ()
{
    BL_PROFILE("incflo::InitData()");

    int restart_flag = 0;
    if(m_restart_file.empty())
    {
        // This tells the AmrMesh class not to iterate when creating the initial
        // grid hierarchy
        // SetIterateToFalse();

        // This tells the Cluster routine to use the new chopping routine which
        // rejects cuts if they don't improve the efficiency
        SetUseNewChop();

        // This is an AmrCore member function which recursively makes new levels
        // with MakeNewLevelFromScratch.
        InitFromScratch(m_cur_time);

        if (m_do_initial_proj) {
            InitialProjection();
        }
        if (m_initial_iterations > 0) {
            InitialIterations();
        }

        // Set m_nstep to 0 before entering time loop
        m_nstep = 0;

        // xxxxx TODO averagedown ???

        if (m_check_int > 0) { WriteCheckPointFile(); }
    }
    else
    {
        restart_flag = 1;
        // Read starting configuration from chk file.
        ReadCheckpointFile();
    }

    // Plot initial distribution
    if((m_plot_int > 0 || m_plot_per_exact > 0 || m_plot_per_approx > 0) && !restart_flag)
    {
        WritePlotFile();
        m_last_plt = 0;
    }
    if(m_KE_int > 0 && !restart_flag)
    {
        amrex::Abort("xxxxx m_KE_int todo");
//        amrex::Print() << "Time, Kinetic Energy: " << m_cur_time << ", " << ComputeKineticEnergy() << std::endl;
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

    if (m_verbose > 0 and ParallelDescriptor::IOProcessor()) {
        printGridSummary(amrex::OutStream(), 0, finest_level);
    }
}

void incflo::Evolve()
{
    BL_PROFILE("incflo::Evolve()");

    bool do_not_evolve = ((m_max_step == 0) || ((m_stop_time >= 0.) && (m_cur_time > m_stop_time)) ||
   					     ((m_stop_time <= 0.) && (m_max_step <= 0))) && !m_steady_state;

    while(!do_not_evolve)
    {
        if (m_verbose > 0)
        {
            amrex::Print() << "\n ============   NEW TIME STEP   ============ \n";
        }

        if (m_regrid_int > 0 and m_nstep > 0 and m_nstep%m_regrid_int == 0)
        {
            if (m_verbose > 0) amrex::Print() << "Regriding...\n";
            regrid(0, m_cur_time);
            if (m_verbose > 0 and ParallelDescriptor::IOProcessor()) {
                printGridSummary(amrex::OutStream(), 0, finest_level);
            }
        }

        // Advance to time t + dt
        Advance();
        m_nstep++;
        m_cur_time += m_dt;

        if (writeNow())
        {
            WritePlotFile();
            m_last_plt = m_nstep;
        }

        if(m_check_int > 0 && (m_nstep % m_check_int == 0))
        {
            WriteCheckPointFile();
            m_last_chk = m_nstep;
        }
        
        if(m_KE_int > 0 && (m_nstep % m_KE_int == 0))
        {
            amrex::Print() << "Time, Kinetic Energy: " << m_cur_time << ", " << ComputeKineticEnergy() << std::endl;
        }

        // Mechanism to terminate incflo normally.
        do_not_evolve = (m_steady_state && SteadyStateReached()) ||
                        ((m_stop_time > 0. && (m_cur_time >= m_stop_time - 1.e-12 * m_dt)) ||
                         (m_max_step >= 0 && m_nstep >= m_max_step));
    }

	// Output at the final time
    if( m_check_int > 0 && m_nstep != m_last_chk) {
        WriteCheckPointFile();
    }
    if( (m_plot_int > 0 || m_plot_per_exact > 0 || m_plot_per_approx > 0)
        && m_nstep != m_last_plt)
    {
        WritePlotFile();
    }
}

// Make a new level from scratch using provided BoxArray and DistributionMapping.
// Only used during initialization.
// overrides the pure virtual function in AmrCore
void incflo::MakeNewLevelFromScratch (int lev, Real time, const BoxArray& new_grids,
                                      const DistributionMapping& new_dmap)
{
    BL_PROFILE("incflo::MakeNewLevelFromScratch()");

    if (m_verbose > 0)
    {
        amrex::Print() << "Making new level " << lev << " from scratch" << std::endl;
        if (m_verbose > 2) {
            amrex::Print() << "with BoxArray " << new_grids << std::endl;
        }
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
                                         m_ntrac, nghost_state(),
                                         m_use_godunov,
                                         m_diff_type==DiffusionType::Implicit,
                                         m_advect_tracer));

    m_t_new[lev] = time;
    m_t_old[lev] = time - 1.e200;

    if (m_restart_file.empty()) {
        prob_init_fluid(lev);
    }
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

    if ( m_plot_int > 0 && (m_nstep % m_plot_int == 0) ) 
        write_now = true;

    else if ( m_plot_per_exact  > 0 && (std::abs(std::remainder(m_cur_time, m_plot_per_exact)) < 1.e-12) ) 
        write_now = true;

    else if (m_plot_per_approx > 0.0)
    {
        // Check to see if we've crossed a m_plot_per_approx interval by comparing
        // the number of intervals that have elapsed for both the current
        // time and the time at the beginning of this timestep.

        int num_per_old = (m_cur_time-m_dt) / m_plot_per_approx;
        int num_per_new = (m_cur_time     ) / m_plot_per_approx;

        // Before using these, however, we must test for the case where we're
        // within machine epsilon of the next interval. In that case, increment
        // the counter, because we have indeed reached the next m_plot_per_approx interval
        // at this point.

        const Real eps = std::numeric_limits<Real>::epsilon() * 10.0 * std::abs(m_cur_time);
        const Real next_plot_time = (num_per_old + 1) * m_plot_per_approx;

        if ((num_per_new == num_per_old) && std::abs(m_cur_time - next_plot_time) <= eps)
        {
            num_per_new += 1;
        }

        // Similarly, we have to account for the case where the old time is within
        // machine epsilon of the beginning of this interval, so that we don't double
        // count that time threshold -- we already plotted at that time on the last timestep.

        if ((num_per_new != num_per_old) && std::abs((m_cur_time - m_dt) - next_plot_time) <= eps)
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

Vector<MultiFab*> incflo::get_divtau_old () noexcept
{
    Vector<MultiFab*> r;
    r.reserve(finest_level+1);
    for (int lev = 0; lev <= finest_level; ++lev) {
        r.push_back(&(m_leveldata[lev]->divtau_o));
    }
    return r;
}

Vector<MultiFab*> incflo::get_divtau_new () noexcept
{
    Vector<MultiFab*> r;
    r.reserve(finest_level+1);
    for (int lev = 0; lev <= finest_level; ++lev) {
        r.push_back(&(m_leveldata[lev]->divtau));
    }
    return r;
}

Vector<MultiFab*> incflo::get_laps_old () noexcept
{
    Vector<MultiFab*> r;
    r.reserve(finest_level+1);
    for (int lev = 0; lev <= finest_level; ++lev) {
        r.push_back(&(m_leveldata[lev]->laps_o));
    }
    return r;
}

Vector<MultiFab*> incflo::get_laps_new () noexcept
{
    Vector<MultiFab*> r;
    r.reserve(finest_level+1);
    for (int lev = 0; lev <= finest_level; ++lev) {
        r.push_back(&(m_leveldata[lev]->laps));
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
    if (m_ntrac > 0) {
        MultiFab::Copy(m_leveldata[lev]->tracer_o,
                       m_leveldata[lev]->tracer, 0, 0, m_ntrac, ng);
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
    if (m_ntrac > 0) {
        MultiFab::Copy(m_leveldata[lev]->tracer,
                       m_leveldata[lev]->tracer_o, 0, 0, m_ntrac, ng);
    }
}
