#include "ocean_waves_test_utils.H"

namespace amr_wind_tests {
namespace utils {

void populate_ocean_waves_params()
{
    // Needed for initial conditions
    {
        amrex::ParmParse pp("incflo");
        pp.add("use_godunov", 1);
        amrex::Vector<std::string> physics{"OceanWaves"};
        pp.addarr("physics", physics);

        pp.add("density", 1.0); // Density
        amrex::Vector<amrex::Real> vel{{20.0, 10.0, 0.0}};
        pp.addarr("velocity", vel);
        amrex::Vector<amrex::Real> grav{{0.0, 0.0, -9.81}};
        pp.addarr("gravity", grav);
    }

    // Adjust computational domain to be more like the ocean waves mesh in the z
    // direction
    {
        amrex::ParmParse pp("amr");
        amrex::Vector<int> ncell{{64, 8, 32}};
        pp.addarr("n_cell", ncell);
    }
    {
        amrex::ParmParse pp("geometry");
        amrex::Vector<amrex::Real> probhi{{4.0, 0.5, 1.0}};
        pp.addarr("prob_hi", probhi);
    }
}

} // namespace utils
} // namespace amr_wind_tests
