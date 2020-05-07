#include "mms_test_utils.H"

namespace amr_wind_tests {
namespace utils {
void populate_mms_params()
{
    amrex::ParmParse pp("MMS");
    const std::string masa_name = "navierstokes_3d_incompressible_homogeneous";
    pp.add("masa_name", masa_name);

    // transport
    {
        amrex::ParmParse pp("transport");
        pp.add("viscosity", 1.0);
    }

    // MMS Forcing
    {
        amrex::ParmParse pp("ICNS");
        pp.add("source_terms", "MMSForcing");
    }

    // Needed for initial conditions
    {
        amrex::ParmParse pp("incflo");
        pp.add("probtype", 0);
        pp.add("use_godunov", 1);
        amrex::Vector<std::string> physics{"MMS"};
        pp.addarr("physics", physics);

        amrex::Vector<amrex::Real> grav{{0.0, 0.0, 0.0}};
        pp.addarr("gravity", grav);
    }

    // Computational domain
    {
        amrex::ParmParse pp("amr");
        amrex::Vector<int> ncell{{16, 16, 16}};
        pp.addarr("n_cell", ncell);
    }
    {
        amrex::ParmParse pp("geometry");
        amrex::Vector<amrex::Real> problo{
            {-3.14159265358979323, -3.14159265358979323, -3.14159265358979323}};
        amrex::Vector<amrex::Real> probhi{
            {3.14159265358979323, 3.14159265358979323, 3.14159265358979323}};

        pp.addarr("prob_lo", problo);
        pp.addarr("prob_hi", probhi);
    }
    {
        amrex::ParmParse pp("time");
        pp.add("stop_time", 2.0);
        pp.add("fixed_dt", 0.05);
    }
}
} // namespace utils
} // namespace amr_wind_tests
