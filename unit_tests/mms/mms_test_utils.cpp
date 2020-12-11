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
        amrex::ParmParse pp2("transport");
        pp2.add("viscosity", 1.0);
    }

    // MMS Forcing
    {
        amrex::ParmParse pp2("ICNS");
        pp2.add("source_terms", "MMSForcing");
    }

    // Needed for initial conditions
    {
        amrex::ParmParse pp2("incflo");
        pp2.add("probtype", 0);
        pp2.add("use_godunov", 1);
        amrex::Vector<std::string> physics{"MMS"};
        pp2.addarr("physics", physics);

        amrex::Vector<amrex::Real> grav{{0.0, 0.0, 0.0}};
        pp2.addarr("gravity", grav);
    }

    // Computational domain
    {
        amrex::ParmParse pp2("amr");
        amrex::Vector<int> ncell{{16, 16, 16}};
        pp2.addarr("n_cell", ncell);
    }
    {
        amrex::ParmParse pp2("geometry");
        amrex::Vector<amrex::Real> problo{
            {-3.14159265358979323, -3.14159265358979323, -3.14159265358979323}};
        amrex::Vector<amrex::Real> probhi{
            {3.14159265358979323, 3.14159265358979323, 3.14159265358979323}};

        pp2.addarr("prob_lo", problo);
        pp2.addarr("prob_hi", probhi);
    }
    {
        amrex::ParmParse pp2("time");
        pp2.add("stop_time", 2.0);
        pp2.add("fixed_dt", 0.05);
    }
}
} // namespace utils
} // namespace amr_wind_tests
