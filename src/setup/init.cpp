#include <AMReX_ParmParse.H>
#include <AMReX_BC_TYPES.H>
#include <AMReX_Box.H>

#include <incflo.H>
#include <boundary_conditions_F.H>
#include <setup_F.H>

#ifdef AMREX_USE_EB
#include <embedded_boundaries_F.H>
#endif

void incflo::ReadParameters()
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
    ReadABLParameters();
    
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
        if (cfl > 0.5)
            amrex::Abort("We currently require cfl <= 0.5 with the current advection scheme");

        // Physics
	pp.queryarr("delp", delp, 0, AMREX_SPACEDIM);
	pp.queryarr("gravity", gravity, 0, AMREX_SPACEDIM);

        pp.query("constant_density", constant_density);
        pp.query("advect_tracer"   , advect_tracer);
        pp.query("use_godunov"     , use_godunov);


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
        AMREX_ALWAYS_ASSERT(ro_0 >= 0.0);

        pp.query("ntrac", ntrac);

        if (ntrac < 1)
            amrex::Abort("We currently require at least one tracer");

        // Scalar diffusion coefficients
        mu_s.resize(ntrac);
        for (int i = 0; i < ntrac; i++) mu_s[i] = 0.;
        pp.queryarr("mu_s", mu_s, 0, ntrac );

        amrex::Print() << "Scalar diffusion coefficients " << std::endl;
        for (int i = 0; i < ntrac; i++)
           amrex::Print() << "Tracer" << i << ":" << mu_s[i] << std::endl;

        AMREX_ALWAYS_ASSERT(mu > 0.0);

        // Get cyclicity, (to pass to Fortran)
        Vector<int> is_cyclic(AMREX_SPACEDIM);
        for(int dir = 0; dir < AMREX_SPACEDIM; dir++)
        {
            is_cyclic[dir] = geom[0].isPeriodic(dir);
        }

        // Loads constants given at runtime `inputs` file into the Fortran module "constant"
        fortran_get_data(is_cyclic.dataPtr(),
                         delp.dataPtr(), gravity.dataPtr(), &ro_0, &mu,
                         &ic_u, &ic_v, &ic_w, &ic_p,
                         &n_0, &ntrac, &tau_0, &papa_reg, &eta_0,
                         fluid_model.c_str(), fluid_model.size());
    } // end prefix incflo

    { // Prefix mac
        ParmParse pp_mac("mac");
        pp_mac.query( "mg_verbose"   , mac_mg_verbose );
        pp_mac.query( "mg_cg_verbose", mac_mg_cg_verbose );
        pp_mac.query( "mg_rtol"      , mac_mg_rtol );
        pp_mac.query( "mg_atol"      , mac_mg_atol );
        pp_mac.query( "mg_maxiter"   , mac_mg_maxiter );
        pp_mac.query( "mg_cg_maxiter", mac_mg_cg_maxiter );
        pp_mac.query( "mg_max_coarsening_level", mac_mg_max_coarsening_level );
    } // end prefix mac
}

void incflo::ReadABLParameters()
{
    ParmParse pp("abl");

    // ABL Physics
    pp.query("ntemperature", ntemperature);
    AMREX_ALWAYS_ASSERT(ntemperature > 0);

    temperature_values.resize(ntemperature);
    temperature_heights.resize(ntemperature);
    for (int i = 0; i < ntemperature; i++) {
        temperature_heights[i] = 650.0;
        temperature_values[i] = 300.0;
    }
    
    pp.queryarr("temperature_heights", temperature_heights,0,ntemperature);
    pp.queryarr("temperature_values", temperature_values,0,ntemperature);
    
    pp.query("use_boussinesq", use_boussinesq);
    pp.query("coriolis_effect", coriolis_effect);
    pp.query("abl_forcing", abl_forcing);
    pp.query("sgs_model", sgs_model);
    
    pp.query("cutoff_height",cutoff_height);
    pp.query("Uperiods",Uperiods);
    pp.query("Vperiods",Vperiods);
    pp.query("deltaU",deltaU);
    pp.query("deltaV",deltaV);
    pp.query("zRefHeight",zRefHeight);
    pp.query("theta_amplitude",theta_amplitude);
    pp.query("abl_forcing_height",abl_forcing_height);
    pp.query("kappa",kappa);
    pp.query("surface_roughness_z0",surface_roughness_z0);
    pp.query("corfac",corfac);
    pp.query("latitude",latitude);

    // set the default to be 1/T0 unless it exists and then it will override
    thermalExpansionCoeff = 1.0/temperature_values[0];
    pp.query("thermalExpansionCoeff",thermalExpansionCoeff);

    pp.query("Smagorinsky_Lilly_SGS_constant",Smagorinsky_Lilly_SGS_constant);
    
}
void incflo::ReadIOParameters()
{
    // Prefix amr
    ParmParse pp("amr");

    pp.query("check_file", check_file);
    pp.query("check_int", check_int);
    pp.query("restart", restart_file);

    pp.query("plot_file", plot_file);
    pp.query("plot_int", plot_int);
    pp.query("plot_per", plot_per);

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
    // Set the BC types on domain boundary
    SetBCTypes();

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

    if(probtype==35) spatially_average_quantities_down(false);// fixme this is only done so that vx_mean_ground is computed for shear stress in bc
    
    // Fill boundaries
    incflo_set_density_bcs(cur_time, density);
    incflo_set_tracer_bcs (cur_time, tracer);
    incflo_set_density_bcs(cur_time, density_o);
    incflo_set_velocity_bcs(cur_time, vel, 0);
    incflo_set_tracer_bcs (cur_time, tracer_o);

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
	Real xlo = geom[0].ProbLo(0);
	Real ylo = geom[0].ProbLo(1);
	Real zlo = geom[0].ProbLo(2);

    for(int lev = 0; lev <= max_level; lev++)
    {
        Box domain(geom[lev].Domain());

        Real dx = geom[lev].CellSize(0);
        Real dy = geom[lev].CellSize(1);
        Real dz = geom[lev].CellSize(2);

        // We deliberately don't tile this loop since we will be looping
        //    over bc's on faces and it makes more sense to do this one grid at a time
        for(MFIter mfi(*density[lev], false); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.validbox();
            const Box& sbx = (*density[lev])[mfi].box();
            init_fluid(sbx.loVect(), sbx.hiVect(),
                       bx.loVect(), bx.hiVect(),
                       domain.loVect(), domain.hiVect(),
                       (*p[lev])[mfi].dataPtr(),
                       (*vel[lev])[mfi].dataPtr(),
                       (*density[lev])[mfi].dataPtr(),
                       (*tracer[lev])[mfi].dataPtr(),
                       (*eta[lev])[mfi].dataPtr(),
                       &dx, &dy, &dz,
                       &xlo, &ylo, &zlo, &probtype);
        }

        if(probtype == 35) {
            init_abl(*density[lev],*vel[lev],*tracer[lev],dx,dy,dz);
        }
        
        // Make sure to set periodic bc's
            vel[lev]->FillBoundary(geom[lev].periodicity());
        density[lev]->FillBoundary(geom[lev].periodicity());
         tracer[lev]->FillBoundary(geom[lev].periodicity());

        MultiFab::Copy(    *vel_o[lev],     *vel[lev], 0, 0,     vel[lev]->nComp(),     vel_o[lev]->nGrow());
        MultiFab::Copy(*density_o[lev], *density[lev], 0, 0, density[lev]->nComp(), density_o[lev]->nGrow());
        MultiFab::Copy(* tracer_o[lev],  *tracer[lev], 0, 0,  tracer[lev]->nComp(),  tracer_o[lev]->nGrow());
    }
}

void incflo::SetBCTypes()
{
    for(int lev = 0; lev <= max_level; lev++)
    {
        Real dx = geom[lev].CellSize(0);
        Real dy = geom[lev].CellSize(1);
        Real dz = geom[lev].CellSize(2);
        Real xlen = geom[lev].ProbHi(0) - geom[lev].ProbLo(0);
        Real ylen = geom[lev].ProbHi(1) - geom[lev].ProbLo(1);
        Real zlen = geom[lev].ProbHi(2) - geom[lev].ProbLo(2);
        Box domain(geom[lev].Domain());

        set_bc_type(bc_ilo[lev]->dataPtr(), bc_ihi[lev]->dataPtr(),
                    bc_jlo[lev]->dataPtr(), bc_jhi[lev]->dataPtr(),
                    bc_klo[lev]->dataPtr(), bc_khi[lev]->dataPtr(),
                    domain.loVect(), domain.hiVect(),
                    &dx, &dy, &dz, &xlen, &ylen, &zlen, &nghost);
    }
}

void incflo::SetBackgroundPressure()
{
    Real xlen = geom[0].ProbHi(0) - geom[0].ProbLo(0);
    Real ylen = geom[0].ProbHi(1) - geom[0].ProbLo(1);
    Real zlen = geom[0].ProbHi(2) - geom[0].ProbLo(2);

    int delp_dir;
    set_delp_dir(&delp_dir);

    IntVect press_per = IntVect(geom[0].isPeriodic(0),
                                geom[0].isPeriodic(1),
                                geom[0].isPeriodic(2));

    // Here we set a separate periodicity flag for p0 because when we use
    // pressure drop (delp) boundary conditions we fill all variables *except* p0
    // periodically
    if(delp_dir > -1)
        press_per[delp_dir] = 0;

    p0_periodicity = Periodicity(press_per);

    for(int lev = 0; lev <= max_level; lev++)
    {
        Real dx = geom[lev].CellSize(0);
        Real dy = geom[lev].CellSize(1);
        Real dz = geom[lev].CellSize(2);

        Box domain(geom[lev].Domain());

        // We deliberately don't tile this loop since we will be looping
        //    over bc's on faces and it makes more sense to do this one grid at a time
        for(MFIter mfi(*density[lev], false); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.validbox();

            set_p0(bx.loVect(), bx.hiVect(),
                   domain.loVect(), domain.hiVect(),
                   BL_TO_FORTRAN_ANYD((*p0[lev])[mfi]),
                   &gp0[0],
                   &dx, &dy, &dz, &xlen, &ylen, &zlen,
                   &delp_dir,
                   bc_ilo[lev]->dataPtr(),
                   bc_ihi[lev]->dataPtr(),
                   bc_jlo[lev]->dataPtr(),
                   bc_jhi[lev]->dataPtr(),
                   bc_klo[lev]->dataPtr(),
                   bc_khi[lev]->dataPtr(),
                   &nghost);
        }
        p0[lev]->FillBoundary(p0_periodicity);
    }

    if (probtype == 11)
    {
       gp0[0] = 0.; gp0[1] = 0.; gp0[2] = 0.;
       for(int lev = 0; lev <= max_level; lev++)
          p0[lev]->setVal(0.);
 
       use_boussinesq = true;
    }
}

//
// Perform initial pressure iterations
//
void incflo::InitialIterations()
{
    BL_PROFILE("incflo::InitialIterations()");

    int initialisation = 1;
    bool explicit_diffusion = (m_diff_type == DiffusionType::Explicit);
    ComputeDt(initialisation, explicit_diffusion);

    if(incflo_verbose)
    {
        amrex::Print() << "Doing initial pressure iterations with dt = " << dt << std::endl;
    }

    // Fill ghost cells
    incflo_set_density_bcs(cur_time, density);
    incflo_set_velocity_bcs(cur_time, vel, 0);
    incflo_set_tracer_bcs(cur_time, tracer);

    // Copy vel into vel_o
    for(int lev = 0; lev <= finest_level; lev++)
    {
        MultiFab::Copy(*vel_o[lev], *vel[lev], 0, 0, vel[lev]->nComp(), vel_o[lev]->nGrow());
    }

    for (int iter = 0; iter < initial_iterations; ++iter)
    {
        if(incflo_verbose) amrex::Print() << "\n In initial_iterations: iter = " << iter << "\n";

 	ApplyPredictor();

        for (int lev = 0; lev <= finest_level; lev++)
        {
            // Replace vel, density, tracer  by the original values
            MultiFab::Copy(*    vel[lev],     *vel_o[lev], 0, 0,     vel[lev]->nComp(),     vel[lev]->nGrow());
            MultiFab::Copy(*density[lev], *density_o[lev], 0, 0, density[lev]->nComp(), density[lev]->nGrow());
            MultiFab::Copy( *tracer[lev],  *tracer_o[lev], 0, 0,  tracer[lev]->nComp(),  tracer[lev]->nGrow());
        }

        // Reset the boundary values (necessary if they are time-dependent)
        incflo_set_velocity_bcs(cur_time, vel    , 0);
        incflo_set_density_bcs (cur_time, density);
        incflo_set_tracer_bcs  (cur_time, tracer );
    }
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
    ApplyProjection(cur_time, dummy_dt);

    // Set nstep (initially -1) to 0, so that subsequent call to ApplyProjection()
    // use the correct decomposition.
    nstep = 0;

    // We set p and gp back to zero (p0 may still be still non-zero)
    for (int lev = 0; lev <= finest_level; lev++)
    {
        p[lev]->setVal(0.0);
        gp[lev]->setVal(0.0);
    }

    if (incflo_verbose)
    {
        amrex::Print() << "After initial projection:" << std::endl;
        PrintMaxValues(time);
    }
}
