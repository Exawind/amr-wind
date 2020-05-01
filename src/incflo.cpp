
#include <incflo.H>

#include "ABL.H"
#include "RefinementCriteria.H"
#include "PDEBase.H"

using namespace amrex;

incflo::incflo ()
    : m_sim(*this)
    , m_time(m_sim.time())
    , m_repo(m_sim.repo())
{
    // NOTE: Geometry on all levels has just been defined in the AmrCore
    // constructor. No valid BoxArray and DistributionMapping have been defined.
    // But the arrays for them have been resized.

    m_time.parse_parameters();
    // Read inputs file using ParmParse
    ReadParameters();

    declare_fields();

    init_field_bcs();

    set_background_pressure();
}

incflo::~incflo ()
{}

void incflo::InitData ()
{
    BL_PROFILE("amr-wind::incflo::InitData()")

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
        amrex::Print() << "Creating mesh... " ;
        InitFromScratch(m_time.current_time());
        amrex::Print() << "done" << std::endl;
        if (ParallelDescriptor::IOProcessor()) {
            amrex::Print() << "Grid summary: " << std::endl;
            printGridSummary(amrex::OutStream(), 0, finest_level);
        }
        for (auto& pp: m_sim.physics())
            pp->post_init_actions();

        icns().initialize();
        for (auto& eqn: scalar_eqns()) eqn->initialize();

        if (m_do_initial_proj) {
            InitialProjection();
        }
        if (m_initial_iterations > 0) {
            InitialIterations();
        }

        if (m_time.write_checkpoint()) { WriteCheckPointFile(); }
    }
    else
    {
        restart_flag = 1;
        // Read starting configuration from chk file.
        ReadCheckpointFile();
        if (ParallelDescriptor::IOProcessor()) {
            amrex::Print() << "Grid summary: " << std::endl;
            printGridSummary(amrex::OutStream(), 0, finest_level);
        }

        icns().initialize();
        for (auto& eqn: scalar_eqns()) eqn->initialize();
    }

    // Plot initial distribution
    if(m_time.write_plot_file() && !restart_flag)
    {
        WritePlotFile();
        m_last_plt = 0;
    }
}

void incflo::Evolve()
{
    BL_PROFILE("amr-wind::incflo::Evolve()")

    if (m_KE_int > 0 && m_restart_file.empty()) {
        amrex::Print() << "\nTime, Kinetic Energy: " << m_time.current_time()
                       << ", " << ComputeKineticEnergy() << std::endl;
    }

    while(m_time.new_timestep())
    {
        amrex::Real time0 = amrex::ParallelDescriptor::second();
        if (m_time.do_regrid())
        {
            amrex::Print() << "Regrid mesh ... " ;
            amrex::Real rstart = amrex::ParallelDescriptor::second();
            regrid(0, m_time.current_time());
            amrex::Real rend = amrex::ParallelDescriptor::second() - rstart;
            amrex::Print() << "time elapsed = " << rend << std::endl;
            if (ParallelDescriptor::IOProcessor()) {
                amrex::Print() << "Grid summary: " << std::endl;
                printGridSummary(amrex::OutStream(), 0, finest_level);
            }
            icns().post_regrid_actions();
            for (auto& eqn: scalar_eqns()) eqn->post_regrid_actions();
        }

        // Advance to time t + dt
        amrex::Real time1 = amrex::ParallelDescriptor::second();
        Advance();
        amrex::Print() << std::endl;
        amrex::Real time2 = amrex::ParallelDescriptor::second();

        if (m_time.write_plot_file())
        {
            WritePlotFile();
            m_last_plt = m_time.time_index();
        }

        if(m_time.write_checkpoint())
        {
            WriteCheckPointFile();
            m_last_chk = m_time.time_index();
        }
        
        if(m_KE_int > 0 && (m_time.time_index() % m_KE_int == 0))
        {
            amrex::Print() << "Time, Kinetic Energy: " << m_time.current_time()
                           << ", " << ComputeKineticEnergy() << std::endl;
        }
        amrex::Real time3 = amrex::ParallelDescriptor::second();

        amrex::Print() << "WallClockTime: " << m_time.time_index()
                       << " Solve: " << std::setprecision(8) << (time2 - time1)
                       << " Misc: " << std::setprecision(6) << (time3 - time2)
                       << " Total: " << std::setprecision(8) << (time3 - time0) << std::endl;
    }
    amrex::Print()
        << "\n==============================================================================\n"
        << std::endl;

    // TODO: Fix last checkpoint/plot output
    // Output at the final time
    if( m_time.write_last_checkpoint()) {
        WriteCheckPointFile();
    }
    if( m_time.write_last_plot_file())
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
    BL_PROFILE("amr-wind::incflo::MakeNewLevelFromScratch()")

    if (m_verbose > 0)
    {
        amrex::Print() << "Making new level " << lev << " from scratch" << std::endl;
        if (m_verbose > 2) {
            amrex::Print() << "with BoxArray " << new_grids << std::endl;
        }
    }

    SetBoxArray(lev, new_grids);
    SetDistributionMap(lev, new_dmap);

    m_repo.make_new_level_from_scratch(lev, time, new_grids, new_dmap);

    if (m_restart_file.empty()) {
        for (auto& pp: m_sim.physics()) {
            pp->initialize_fields(lev, Geom(lev));
        }
    }
}
