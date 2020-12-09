#include "aw_test_utils/AmrexTestEnv.H"
#include "pp_utils.H"
#include "AMReX_ParmParse.H"

extern amr_wind_tests::AmrexTestEnv* utest_env;

namespace amr_wind_tests {
namespace pp_utils {

bool has_managed_memory()
{
#if defined(AMREX_USE_HIP)
    return false;
#else
    return utest_env->has_managed_memory();
#endif
}

void default_time_inputs()
{
    amrex::ParmParse pp("time");
    pp.add("stop_time", 2.0);
    pp.add("max_step", 10);
    pp.add("fixed_dt", 0.1);
    pp.add("init_shrink", 0.1);
    pp.add("cfl", 0.5);
    pp.add("verbose", -1);
    pp.add("regrid_interval", 3);
    pp.add("plot_interval", 1);
    pp.add("checkpoint_interval", 2);
}

void default_mesh_inputs()
{
    {
        amrex::ParmParse pp("amr");
        amrex::Vector<int> ncell{{8, 8, 8}};

        pp.add("verbose", 0);
        pp.addarr("n_cell", ncell);
        pp.add("max_level", 0);
    }

    {
        amrex::ParmParse pp("geometry");
        amrex::Vector<amrex::Real> problo{{0.0, 0.0, 0.0}};
        amrex::Vector<amrex::Real> probhi{{8.0, 8.0, 8.0}};
        amrex::Vector<int> periodic{{1, 1, 1}};

        pp.addarr("prob_lo", problo);
        pp.addarr("prob_hi", probhi);
        pp.addarr("is_periodic", periodic);
    }

    {
        amrex::ParmParse pp("incflo");
        pp.add("probtype", 0);
    }
}

} // namespace pp_utils
} // namespace amr_wind_tests
