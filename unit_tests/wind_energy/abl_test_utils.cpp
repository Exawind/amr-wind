#include "abl_test_utils.H"

namespace amr_wind_tests {
namespace utils {

void populate_abl_params()
{
    amrex::ParmParse pp("ABL");

    // Initial conditions (Temperature)
    amrex::Vector<amrex::Real> theights{{0.0, 650.0, 750.0, 1000.0}};
    amrex::Vector<amrex::Real> tvalues{{300.0, 300.0, 308.0, 308.75}};
    pp.addarr("temperature_heights", theights);
    pp.addarr("temperature_values", tvalues);
    pp.add("perturb_ref_height", 50.0);
    pp.add("reference_temperature", 300.0);

    // Body force
    {
        amrex::ParmParse pp2("BodyForce");
        pp2.add("type", std::string("oscillatory"));
        amrex::Vector<amrex::Real> source_mag{{1.0, 2.0, 3.0}};
        pp2.addarr("magnitude", source_mag);
        pp2.add("angular_frequency", 1.0);
    }

    // Boussinesq Buoyancy
    {
        amrex::ParmParse pp2("BoussinesqBuoyancy");
        pp2.add("reference_temperature", 300.0);
    }

    // ABL Forcing
    {
        amrex::ParmParse pp2("ABLForcing");
        pp2.add("abl_forcing_height", 90.0);
    }

    // Geostrophic Forcing
    {
        amrex::ParmParse pp2("GeostrophicForcing");
        amrex::Vector<amrex::Real> gwind{{10.0, 6.0, 0.0}};
        pp2.addarr("geostrophic_wind", gwind);
    }

    // Coriolis term
    {
        amrex::ParmParse pp2("CoriolisForcing");
        pp2.add("latitude", 45.0);
    }

    pp.add("kappa", 0.41);
    pp.add("surface_roughness_z0", 0.1);

    // Needed for initial conditions
    {
        amrex::ParmParse pp2("incflo");
        pp2.add("probtype", 0);
        pp2.add("use_godunov", 1);
        amrex::Vector<std::string> physics{"ABL"};
        pp2.addarr("physics", physics);

        pp2.add("density", 1.0); // Density
        amrex::Vector<amrex::Real> vel{{20.0, 10.0, 0.0}};
        pp2.addarr("velocity", vel);
        amrex::Vector<amrex::Real> grav{{0.0, 0.0, -9.81}};
        pp2.addarr("gravity", grav);
    }

    // Adjust computational domain to be more like ABL mesh in the z direction
    {
        amrex::ParmParse pp2("amr");
        amrex::Vector<int> ncell{{8, 8, 64}};
        pp2.addarr("n_cell", ncell);
    }
    {
        amrex::ParmParse pp2("geometry");
        amrex::Vector<amrex::Real> probhi{{120.0, 120.0, 1000.0}};

        pp2.addarr("prob_hi", probhi);
    }
}

} // namespace utils
} // namespace amr_wind_tests
