#include "pp_utils.H"
#include "AMReX_ParmParse.H"

namespace amr_wind_tests {
namespace pp_utils {

void default_time_inputs()
{
    {
        amrex::ParmParse pp;
        pp.add("verbose", -1);
        pp.add("stop_time", 2.0);
        pp.add("max_step", 10.0);
    }

    {
        amrex::ParmParse pp("amr");
        pp.add("regrid_int", 3);
        pp.add("plot_int", 1);
        pp.add("check_int", 2);
    }

    {
        amrex::ParmParse pp("incflo");
        pp.add("fixed_dt", 0.1);
        pp.add("cfl", 0.5);
        pp.add("m_init_shrink", 0.1);
    }
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
        amrex::Vector<amrex::Real> problo {{0.0, 0.0, 0.0}};
        amrex::Vector<amrex::Real> probhi {{8.0, 8.0, 8.0}};
        amrex::Vector<int> periodic{{1, 1, 1}};

        pp.addarr("prob_lo", problo);
        pp.addarr("prob_hi", probhi);
        pp.addarr("is_periodic", periodic);
    }
}

}
}
