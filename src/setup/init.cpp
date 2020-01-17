// #include <AMReX_ParmParse.H>
#include <AMReX_BC_TYPES.H>
#include <incflo.H>

using namespace amrex;

void incflo::ReadParameters ()
{
    {
        // Variables without prefix in inputs file
	ParmParse pp;

	pp.query("stop_time", stop_time);
	pp.query("max_step", max_step);
	pp.query("steady_state", steady_state);
    }

    ReadIOParameters();
    ReadRheologyParameters();

    { // Prefix amr
 	ParmParse pp("amr");

	pp.query("regrid_int", regrid_int);
#ifdef AMREX_USE_EB
        pp.query("refine_cutcells", refine_cutcells);
#endif

        pp.query("KE_int", KE_int);

    } // end prefix amr

    { // Prefix incflo
	ParmParse pp("incflo");

        pp.query("verbose", incflo_verbose);

	pp.query("steady_state_tol", steady_state_tol);
        pp.query("initial_iterations", initial_iterations);
        pp.query("do_initial_proj", do_initial_proj);

	pp.query("fixed_dt", fixed_dt);
	pp.query("cfl", cfl);

        // This will multiply the time-step in the very first step only
	pp.query("init_shrink", init_shrink);
        if (init_shrink > 1.0)
            amrex::Abort("We require init_shrink <= 1.0");

        // Physics
	pp.queryarr("delp", delp, 0, AMREX_SPACEDIM);
	pp.queryarr("gravity", gravity, 0, AMREX_SPACEDIM);

        pp.query("constant_density", constant_density);
        pp.query("advect_tracer"   , advect_tracer);
        pp.query("test_tracer_conservation" , test_tracer_conservation);
        pp.query("use_godunov"        , use_godunov);
        pp.query("use_forces_in_trans", use_forces_in_trans);

        // The default for diffusion_type is 2, i.e. the default m_diff_type is DiffusionType::Implicit
        pp.query("diffusion_type", diffusion_type);
        if (diffusion_type == 0)
           m_diff_type = DiffusionType::Explicit;
        else if (diffusion_type == 1)
           m_diff_type = DiffusionType::Crank_Nicolson;
        else if (diffusion_type == 2)
           m_diff_type = DiffusionType::Implicit;
        else
            amrex::Abort("We currently require diffusion_type = 0 for explicit, 1 for Crank-Nicolson or 2 for implicit");

        if (!use_godunov && cfl > 0.5)
            amrex::Abort("We currently require cfl <= 0.5 when using the MOL advection scheme");
        if (use_godunov && cfl > 1.0)
            amrex::Abort("We currently require cfl <= 1.0 when using the Godunov advection scheme");

        AMREX_ALWAYS_ASSERT(ro_0 >= 0.0);

        // Initial conditions
        pp.query("probtype", probtype);
        pp.query("ic_u", ic_u);
        pp.query("ic_v", ic_v);
        pp.query("ic_w", ic_w);
        pp.query("ic_p", ic_p);

        // Viscosity (if constant)
        pp.query("mu", mu);

        // Density (if constant)
        pp.query("ro_0", ro_0);

        pp.query("ntrac", ntrac);

        if (ntrac <= 0) advect_tracer = 0;

        if (ntrac < 1) {
            amrex::Abort("We currently require at least one tracer");
        }

        // Scalar diffusion coefficients
        mu_s.resize(ntrac);
        for (int i = 0; i < ntrac; i++) mu_s[i] = 0.;
        pp.queryarr("mu_s", mu_s, 0, ntrac );

        amrex::Print() << "Scalar diffusion coefficients " << std::endl;
        for (int i = 0; i < ntrac; i++) {
            amrex::Print() << "Tracer" << i << ":" << mu_s[i] << std::endl;
        }

        // AMREX_ALWAYS_ASSERT(mu > 0.0);
    } // end prefix incflo

    { // Prefix mac
        ParmParse pp_mac("mac_proj");
        pp_mac.query( "mg_verbose"   , mac_mg_verbose );
        pp_mac.query( "mg_cg_verbose", mac_mg_cg_verbose );
        pp_mac.query( "mg_rtol"      , mac_mg_rtol );
        pp_mac.query( "mg_atol"      , mac_mg_atol );
        pp_mac.query( "mg_maxiter"   , mac_mg_maxiter );
        pp_mac.query( "mg_cg_maxiter", mac_mg_cg_maxiter );
        pp_mac.query( "mg_max_coarsening_level", mac_mg_max_coarsening_level );
    } // end prefix mac
}

void incflo::ReadIOParameters()
{
    // Prefix amr
    ParmParse pp("amr");

    pp.query("check_file", check_file);
    pp.query("check_int", check_int);
    pp.query("restart", restart_file);

    pp.query("plot_file", plot_file);
    pp.query("plot_int"       , plot_int);
    pp.query("plot_per_exact" , plot_per_exact);
    pp.query("plot_per_approx", plot_per_approx);

    if ( (plot_int       > 0 && plot_per_exact  > 0) ||
         (plot_int       > 0 && plot_per_approx > 0) ||
         (plot_per_exact > 0 && plot_per_approx > 0) )
       amrex::Abort("Must choose only one of plot_int or plot_per_exact or plot_per_approx");

    // The plt_ccse_regtest resets the defaults,
    //     but we can over-ride those below
    int plt_ccse_regtest = 0;
    pp.query("plt_ccse_regtest", plt_ccse_regtest);

    if (plt_ccse_regtest != 0)
    {
        plt_velx       = 1;
        plt_vely       = 1;
        plt_velz       = 1;
        plt_gpx        = 1;
        plt_gpy        = 1;
        plt_gpz        = 1;
        plt_rho        = 1;
        plt_tracer     = 1;
        plt_p          = 0;
        plt_eta        = 0;
        plt_vort       = 0;
        plt_strainrate = 0;
        plt_stress     = 0;
        plt_divu       = 0;
        plt_vfrac      = 0;
    }

    // Which variables to write to plotfile
    pltVarCount = 0;

    pp.query("plt_velx",       plt_velx  );
    pp.query("plt_vely",       plt_vely  );
    pp.query("plt_velz",       plt_velz  );

    pp.query("plt_gpx",        plt_gpx );
    pp.query("plt_gpy",        plt_gpy );
    pp.query("plt_gpz",        plt_gpz );

    pp.query("plt_rho",        plt_rho   );
    pp.query("plt_tracer",     plt_tracer);
    pp.query("plt_p",          plt_p     );
    pp.query("plt_eta",        plt_eta   );
    pp.query("plt_vort",       plt_vort  );
    pp.query("plt_strainrate", plt_strainrate);
    pp.query("plt_stress"    , plt_stress);
    pp.query("plt_divu",       plt_divu  );
    pp.query("plt_vfrac",      plt_vfrac );

    // Special test for CCSE regression test. Override all individual
    // flags and save all data to plot file.

    // Count the number of variables to save.
    if(plt_velx       == 1) pltVarCount += 1;
    if(plt_vely       == 1) pltVarCount += 1;
    if(plt_velz       == 1) pltVarCount += 1;
    if(plt_gpx        == 1) pltVarCount += 1;
    if(plt_gpy        == 1) pltVarCount += 1;
    if(plt_gpz        == 1) pltVarCount += 1;
    if(plt_rho        == 1) pltVarCount += 1;
    if(plt_tracer     == 1) pltVarCount += ntrac;
    if(plt_p          == 1) pltVarCount += 1;
    if(plt_eta        == 1) pltVarCount += 1;
    if(plt_vort       == 1) pltVarCount += 1;
    if(plt_strainrate == 1) pltVarCount += 1;
    if(plt_stress     == 1) pltVarCount += 1;
    if(plt_divu       == 1) pltVarCount += 1;
    if(plt_vfrac      == 1) pltVarCount += 1;
}

void incflo::PostInit(int restart_flag)
{
    // Init nodal and diffusion solvers (for now only diffusion)
    // (Note we must do this *after* setting the bc types above)
    incflo_init_solvers();

    // Initial fluid arrays: pressure, velocity, density, viscosity
    if(!restart_flag)
    {
        InitFluid();
    }

    // Set the background pressure and gradients in "DELP" cases
    SetBackgroundPressure();

    // Fill boundaries
    incflo_set_density_bcs (cur_time, density);
    incflo_set_tracer_bcs  (cur_time, tracer);
    incflo_set_density_bcs (cur_time, density_o);
    incflo_set_tracer_bcs  (cur_time, tracer_o);
    incflo_set_velocity_bcs(cur_time, vel);

    setup_level_mask();

    // Project the initial velocity field to make it divergence free
    // Perform initial iterations to find pressure distribution
    if(!restart_flag)
    {
        if (do_initial_proj)
            InitialProjection();
        if (initial_iterations > 0)
            InitialIterations();
    }
}

void incflo::setup_level_mask(){

     BL_PROFILE("incflo::setup_level_mask");

     for(int lev=0;lev<finest_level;++lev) {
        *level_mask[lev] = makeFineMask(grids[lev],dmap[lev], grids[lev+1], IntVect(2), 1, 0);
    }
}

void incflo::InitFluid()
{
}

void incflo::SetBCTypes()
{
}

void incflo::SetBackgroundPressure()
{
}

//
// Perform initial pressure iterations
//
void incflo::InitialIterations ()
{
    BL_PROFILE("incflo::InitialIterations()");

    int initialisation = 1;
    bool explicit_diffusion = (m_diff_type == DiffusionType::Explicit);
    ComputeDt(initialisation, explicit_diffusion);

    if(incflo_verbose)
    {
        amrex::Print() << "Doing initial pressure iterations with dt = " << dt << std::endl;
    }

    copy_from_new_to_old_velocity();
    copy_from_new_to_old_density();
    copy_from_new_to_old_tracer();
    for(int lev = 0; lev <= finest_level; ++lev) t_old[lev] = t_new[lev];

    int ng = nghost_state();
    for (int lev = 0; lev <= finest_level; ++lev) {
        fillpatch_velocity(lev, t_old[lev], m_leveldata[lev]->velocity_o, ng);
        fillpatch_density(lev, t_old[lev], m_leveldata[lev]->density_o, ng);
    }

    // No need to advect tracer in initial pressure iterations
    auto advect_tracer_save = advect_tracer;
    advect_tracer = false;

    for (int iter = 0; iter < initial_iterations; ++iter)
    {
        if(incflo_verbose) amrex::Print() << "\n In initial_iterations: iter = " << iter << "\n";

 	ApplyPredictor(true);

        copy_from_old_to_new_velocity();
        copy_from_old_to_new_density();

        if (use_godunov) {
    amrex::VisMF::Write(m_leveldata[0]->velocity, "vel");
    amrex::VisMF::Write(m_leveldata[0]->density, "rho");
    amrex::VisMF::Write(m_leveldata[0]->tracer, "tra");
    amrex::VisMF::Write(m_leveldata[0]->p, "p");

        amrex::Abort("xxxxx after first pressure iteration");
    }

    }

    advect_tracer = advect_tracer_save;

    // Set nstep to 0 before entering time loop
    nstep = 0;
}

// Project velocity field to make sure initial velocity is divergence-free
void incflo::InitialProjection()
{
    BL_PROFILE("incflo::InitialProjection()");

    Real time = 0.0;

    if (incflo_verbose)
    {
        amrex::Print() << "Initial projection:" << std::endl;
        PrintMaxValues(time);
    }

    Real dummy_dt = 1.0;
    bool incremental = false;
    ApplyProjection(cur_time, dummy_dt, incremental);

    // We set p and gp back to zero (p0 may still be still non-zero)
    for (int lev = 0; lev <= finest_level; lev++)
    {
        m_leveldata[lev]->p.setVal(0.0);
        m_leveldata[lev]->gp.setVal(0.0);
    }

    if (incflo_verbose)
    {
        amrex::Print() << "After initial projection:" << std::endl;
        PrintMaxValues(time);
    }
}
