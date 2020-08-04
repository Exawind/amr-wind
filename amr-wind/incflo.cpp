
#include "amr-wind/incflo.H"

#include "amr-wind/wind_energy/ABL.H"
#include "amr-wind/utilities/tagging/RefinementCriteria.H"
#include "amr-wind/equation_systems/PDEBase.H"
#include "amr-wind/utilities/IOManager.H"
#include "amr-wind/utilities/PostProcessing.H"
#include "amr-wind/overset/OversetManager.H"

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

    init_physics_and_pde();

    set_background_pressure();
}

incflo::~incflo ()
{}

/** Initialize AMR mesh data structure.
 *
 *  Calls the AmrMesh/AmrCore functions to initialize the mesh. For restarted
 *  simulations, it loads the checkpoint file.
 */
void incflo::init_mesh()
{
    BL_PROFILE("amr-wind::incflo::init_mesh");
    // Initialize I/O manager to enable restart and outputs
    auto& io_mgr = m_sim.io_manager();
    io_mgr.initialize_io();

    if (!io_mgr.is_restart()) {
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
    } else {
        // Read starting configuration from chk file.
        ReadCheckpointFile();
        if (ParallelDescriptor::IOProcessor()) {
            amrex::Print() << "Grid summary: " << std::endl;
            printGridSummary(amrex::OutStream(), 0, finest_level);
        }
    }
}

/** Initialize AMR-Wind data structures after mesh has been created.
 *
 *  Modules initialized:
 *    - Registered \ref physics classes
 *    - Registered \ref PDE systems
 *    - Registered post-processing classes
 */
void incflo::init_amr_wind_modules()
{
    BL_PROFILE("amr-wind::incflo::init_amr_wind_modules");
    for (auto& pp: m_sim.physics())
        pp->post_init_actions();

    icns().initialize();
    for (auto& eqn: scalar_eqns()) eqn->initialize();

    if (m_sim.has_overset()) m_sim.overset_manager()->post_init_actions();
    m_sim.post_manager().initialize();
}

/** Initialize flow-field before performing time-integration.
 *
 *  This method calls the incflo::InitialProjection step to ensure that the
 *  velocity field is divergence free. Then it performs a user-defined number of
 *  initial iterations to compute \f$p^{n - 1/2}\f$ necessary for advancing into
 *  the first timestep. These actions are only performed when starting a new
 *  simulation, but not for a restarted simulation.
 */
void incflo::prepare_for_time_integration()
{
    BL_PROFILE("amr-wind::incflo::prepare_for_time_integration");
    // Don't perform initial work if this is a restart
    if (m_sim.io_manager().is_restart()) return;

    if (m_do_initial_proj) {
        InitialProjection();
    }
    if (m_initial_iterations > 0) {
        InitialIterations();
    }

    // Plot initial distribution
    if(m_time.write_plot_file())
        m_sim.io_manager().write_plot_file();
    if (m_time.write_checkpoint())
        m_sim.io_manager().write_checkpoint_file();
}

/** Perform all initialization actions for AMR-Wind.
 *
 *  This is a wrapper method that calls the following methods to perform the actual work.
 *
 *  \callgraph
 */
void incflo::InitData ()
{
    BL_PROFILE("amr-wind::incflo::InitData()");

    init_mesh();
    init_amr_wind_modules();
    prepare_for_time_integration();
}

/** Perform regrid actions at a given timestep.
 *
 *  \return Flag indicating if regrid was performed
 */
bool incflo::regrid_and_update()
{
    BL_PROFILE("amr-wind::incflo::regrid_and_update");

    if (m_time.do_regrid()) {
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
        for (auto& pp : m_sim.physics()) pp->post_regrid_actions();

        if (m_sim.has_overset())
            m_sim.overset_manager()->post_regrid_actions();
    }

    return m_time.do_regrid();
}

/** Perform actions after a timestep
 *
 *  Performs any necessary I/O before advancing to next timestep
 */
void incflo::post_advance_work()
{
    BL_PROFILE("amr-wind::incflo::post_advance_work");
    if (m_time.write_plot_file()) {
        m_sim.io_manager().write_plot_file();
    }

    if (m_time.write_checkpoint()) {
        m_sim.io_manager().write_checkpoint_file();
    }

    if (m_KE_int > 0 && (m_time.time_index() % m_KE_int == 0)) {
        amrex::Print() << "Time, Kinetic Energy: " << m_time.new_time() << ", "
                       << ComputeKineticEnergy() << std::endl;
    }
}

/** Perform time-integration for user-defined time or timesteps.
 *
 *  \callgraph
 */
void incflo::Evolve()
{
    BL_PROFILE("amr-wind::incflo::Evolve()");

    if (m_KE_int > 0 && m_sim.io_manager().restart_file().empty()) {
        amrex::Print() << "\nTime, Kinetic Energy: " << m_time.new_time()
                       << ", " << ComputeKineticEnergy() << std::endl;
    }

    while(m_time.new_timestep()) {
        amrex::Real time0 = amrex::ParallelDescriptor::second();
        regrid_and_update();

        // Advance to time t + dt
        amrex::Real time1 = amrex::ParallelDescriptor::second();
        advance();
        amrex::Print() << std::endl;
        amrex::Real time2 = amrex::ParallelDescriptor::second();
        post_advance_work();
        amrex::Real time3 = amrex::ParallelDescriptor::second();

        amrex::Print() << "WallClockTime: " << m_time.time_index()
                       << " Solve: " << std::setprecision(8) << (time2 - time1)
                       << " Misc: " << std::setprecision(6) << (time3 - time2)
                       << " Total: " << std::setprecision(8) << (time3 - time0) << std::endl;
    }
    amrex::Print()
        << "\n==============================================================================\n"
        << std::endl;

    // Output at final time
    if( m_time.write_last_plot_file())
        m_sim.io_manager().write_plot_file();

    if( m_time.write_last_checkpoint())
        m_sim.io_manager().write_checkpoint_file();

}

// Make a new level from scratch using provided BoxArray and DistributionMapping.
// Only used during initialization.
// overrides the pure virtual function in AmrCore
void incflo::MakeNewLevelFromScratch (int lev, Real time, const BoxArray& new_grids,
                                      const DistributionMapping& new_dmap)
{
    BL_PROFILE("amr-wind::incflo::MakeNewLevelFromScratch()");

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

    for (auto& pp: m_sim.physics()) {
        pp->initialize_fields(lev, Geom(lev));
    }
}

void incflo::init_physics_and_pde()
{
    {
        // Query and activate overset before initializing PDEs and physics
        amrex::ParmParse pp("incflo");
        bool activate_overset = false;
        pp.query("activate_overset", activate_overset);

        if (activate_overset) m_sim.activate_overset();
    }

    auto& pde_mgr = m_sim.pde_manager();

    // Always register incompressible Navier-Stokes equation
    pde_mgr.register_icns();

    // Register density first so that we can compute its `n+1/2` state before
    // other scalars attempt to use it in their computations.
    if (!m_constant_density) {
        if (pde_mgr.scalar_eqns().size() > 0)
            amrex::Abort(
                "For non-constant density, it must be the first equation "
                "registered for the scalar equations");
        pde_mgr.register_transport_pde("Density");
    }

    m_sim.init_physics();
    m_sim.create_turbulence_model();
}
