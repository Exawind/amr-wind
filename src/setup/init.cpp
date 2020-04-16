// #include <AMReX_ParmParse.H>
#include <AMReX_BC_TYPES.H>
#include <incflo.H>
#include <cmath>


#include "Physics.H"
#include "ABL.H"
#include "BoussinesqBubble.H"
#include "RefinementCriteria.H"
#include "CartBoxRefinement.H"

using namespace amrex;

void incflo::ReadParameters ()
{
    ReadIOParameters();
    ReadRheologyParameters();

    { // Prefix amr
        ParmParse pp("amr");
        pp.query("KE_int", m_KE_int);
    } // end prefix amr

    { // Prefix incflo
        ParmParse pp("incflo");

        pp.query("verbose", m_verbose);

        pp.query("initial_iterations", m_initial_iterations);
        pp.query("do_initial_proj", m_do_initial_proj);

        // Physics
        pp.queryarr("delp", m_delp, 0, AMREX_SPACEDIM);
        pp.queryarr("gravity", m_gravity, 0, AMREX_SPACEDIM);

        pp.query("constant_density"         , m_constant_density);
        pp.query("advect_tracer"            , m_advect_tracer);
        pp.query("test_tracer_conservation" , m_test_tracer_conservation);

        // Godunov-related flags
        pp.query("use_godunov"                      , m_use_godunov);
        pp.query("use_ppm"                          , m_godunov_ppm);
        pp.query("godunov_use_forces_in_trans"      , m_godunov_use_forces_in_trans);

        // The default for diffusion_type is 2, i.e. the default m_diff_type is DiffusionType::Implicit
        int diffusion_type = 2;
        pp.query("diffusion_type", diffusion_type);
        if (diffusion_type == 0) {
            m_diff_type = DiffusionType::Explicit;
        } else if (diffusion_type == 1) {
            m_diff_type = DiffusionType::Crank_Nicolson;
        } else if (diffusion_type == 2) {
            m_diff_type = DiffusionType::Implicit;
        } else {
            amrex::Abort("We currently require diffusion_type = 0 for explicit, 1 for Crank-Nicolson or 2 for implicit");
        }

        if (!m_use_godunov && m_time.max_cfl() > 0.5) {
            amrex::Abort("We currently require cfl <= 0.5 when using the MOL advection scheme");
        }
        if (m_use_godunov && m_time.max_cfl() > 1.0) {
            amrex::Abort("We currently require cfl <= 1.0 when using the Godunov advection scheme");
        }

        // Initial conditions
        pp.query("probtype", m_probtype);
        pp.query("ic_u", m_ic_u);
        pp.query("ic_v", m_ic_v);
        pp.query("ic_w", m_ic_w);
        pp.query("ic_p", m_ic_p);

        // Viscosity (if constant)
        pp.query("mu", m_mu);

        // Density (if constant)
        pp.query("ro_0", m_ro_0);
        AMREX_ALWAYS_ASSERT(m_ro_0 >= 0.0);

        pp.query("ntrac", m_ntrac);

        if (m_ntrac <= 0) m_advect_tracer = 0;

        if (m_ntrac < 1) {
            amrex::Abort("We currently require at least one tracer");
        }

        // Scalar diffusion coefficients
        m_mu_s.resize(m_ntrac, 0.0);
        pp.queryarr("mu_s", m_mu_s, 0, m_ntrac );

        amrex::Print() << "Scalar diffusion coefficients " << std::endl;
        for (int i = 0; i < m_ntrac; i++) {
            amrex::Print() << "Tracer" << i << ":" << m_mu_s[i] << std::endl;
        }
    } // end prefix incflo

    // FIXME: clean up WIP logic
    if (m_probtype == 35) {
        m_physics.emplace_back(new amr_wind::ABL(m_time, this));
    }
    
    if (m_probtype == 11) {
        m_physics.emplace_back(new amr_wind::BoussinesqBubble(this));
    }

    {
        // tagging options
        ParmParse pp("tagging");
        bool static_refine = false;
        pp.query("static_refinement", static_refine);
        if (static_refine) {
            std::unique_ptr<amr_wind::CartBoxRefinement> obj(new amr_wind::CartBoxRefinement);
            obj->initialize(*this);

            m_refine_criteria.push_back(std::move(obj));
        }
    }

}

void incflo::ReadIOParameters()
{
    // Prefix amr
    ParmParse pp("amr");

    pp.query("check_file", m_check_file);
    pp.query("restart", m_restart_file);

    pp.query("plot_file", m_plot_file);

    // The plt_ccse_regtest resets the defaults,
    //     but we can over-ride those below
    int plt_ccse_regtest = 0;
    pp.query("plt_ccse_regtest", plt_ccse_regtest);

    if (plt_ccse_regtest != 0)
    {
        m_plt_velx       = 1;
        m_plt_vely       = 1;
        m_plt_velz       = 1;
        m_plt_gpx        = 1;
        m_plt_gpy        = 1;
        m_plt_gpz        = 1;
        m_plt_rho        = 1;
        m_plt_tracer     = 1;
        m_plt_p          = 0;
        m_plt_eta        = 0;
        m_plt_vort       = 0;
        m_plt_strainrate = 0;
        m_plt_stress     = 0;
        m_plt_divu       = 0;
        m_plt_vfrac      = 0;
    }

    // Which variables to write to plotfile

    pp.query("plt_velx",       m_plt_velx  );
    pp.query("plt_vely",       m_plt_vely  );
    pp.query("plt_velz",       m_plt_velz  );

    pp.query("plt_gpx",        m_plt_gpx );
    pp.query("plt_gpy",        m_plt_gpy );
    pp.query("plt_gpz",        m_plt_gpz );

    pp.query("plt_rho",        m_plt_rho   );
    pp.query("plt_tracer",     m_plt_tracer);
    pp.query("plt_p",          m_plt_p     );
    pp.query("plt_eta",        m_plt_eta   );
    pp.query("plt_vort",       m_plt_vort  );
    pp.query("plt_strainrate", m_plt_strainrate);
    pp.query("plt_stress"    , m_plt_stress);
    pp.query("plt_divu",       m_plt_divu  );
    pp.query("plt_vfrac",      m_plt_vfrac );

    pp.query("plt_forcing",    m_plt_forcing );
}

//
// Perform initial pressure iterations
//
void incflo::InitialIterations ()
{
    BL_PROFILE("incflo::InitialIterations()");

    bool explicit_diffusion = (m_diff_type == DiffusionType::Explicit);
    ComputeDt(explicit_diffusion);

    if (m_verbose)
    {
        amrex::Print() << "Doing initial pressure iterations with dt = "
                       << m_time.deltaT() << std::endl;
    }

    auto& vel = velocity();
    auto& rho = density();
    auto& trac = tracer();

    vel.copy_state(amr_wind::FieldState::Old, amr_wind::FieldState::New);
    rho.copy_state(amr_wind::FieldState::Old, amr_wind::FieldState::New);
    trac.copy_state(amr_wind::FieldState::Old, amr_wind::FieldState::New);

    vel.state(amr_wind::FieldState::Old).fillpatch(m_time.current_time());
    rho.state(amr_wind::FieldState::Old).fillpatch(m_time.current_time());
    trac.state(amr_wind::FieldState::Old).fillpatch(m_time.current_time());

    for (int iter = 0; iter < m_initial_iterations; ++iter)
    {
        if (m_verbose) amrex::Print() << "In initial_iterations: iter = " << iter << "\n";

        // fixme turn this on later and delete stuff in ABL.cpp but will have to rebless gold files
//        for (auto& pp: m_physics)
//            pp->pre_advance_work();

        ApplyPredictor(true);

        vel.copy_state(amr_wind::FieldState::New, amr_wind::FieldState::Old);
        rho.copy_state(amr_wind::FieldState::New, amr_wind::FieldState::Old);
        trac.copy_state(amr_wind::FieldState::New, amr_wind::FieldState::Old);
    }
}

// Project velocity field to make sure initial velocity is divergence-free
void incflo::InitialProjection()
{
    BL_PROFILE("incflo::InitialProjection()");

    Real time = 0.0;

    if (m_verbose)
    {
        amrex::Print() << "Initial projection:" << std::endl;
        PrintMaxValues(time);
    }

    Real dummy_dt = 1.0;
    bool incremental = false;
    ApplyProjection(m_repo.get_field("density", amr_wind::FieldState::New).vec_const_ptrs(),
                    m_time.current_time(), dummy_dt, incremental);

    // We set p and gp back to zero (p0 may still be still non-zero)
    pressure().setVal(0.0);
    grad_p().setVal(0.0);

    if (m_verbose)
    {
        amrex::Print() << "After initial projection:" << std::endl;
        PrintMaxValues(time);
    }
}
