#include "abl_test_utils.H"

namespace amr_wind_tests::utils {

void populate_abl_params()
{
    // Initial conditions (Temperature)
    {
        amrex::ParmParse pp("ABL");
        amrex::Vector<amrex::Real> theights{{0.0, 650.0, 750.0, 1000.0}};
        amrex::Vector<amrex::Real> tvalues{{300.0, 300.0, 308.0, 308.75}};
        pp.addarr("temperature_heights", theights);
        pp.addarr("temperature_values", tvalues);
        pp.add("perturb_ref_height", 50.0);
        pp.add("reference_temperature", 300.0);
        pp.add("kappa", 0.41);
        pp.add("surface_roughness_z0", 0.1);
    }

    // Body force
    {
        amrex::ParmParse pp("BodyForce");
        pp.add("type", std::string("oscillatory"));
        amrex::Vector<amrex::Real> source_mag{{1.0, 2.0, 3.0}};
        pp.addarr("magnitude", source_mag);
        pp.add("angular_frequency", 1.0);
    }

    // Boussinesq Buoyancy
    {
        amrex::ParmParse pp("BoussinesqBuoyancy");
        pp.add("reference_temperature", 300.0);
    }

    // ABL Forcing
    {
        amrex::ParmParse pp("ABLForcing");
        pp.add("abl_forcing_height", 90.0);
    }

    // Geostrophic Forcing
    {
        amrex::ParmParse pp("GeostrophicForcing");
        amrex::Vector<amrex::Real> gwind{{10.0, 6.0, 0.0}};
        pp.addarr("geostrophic_wind", gwind);
        pp.add("latitude", 45.0);
    }

    // Coriolis term
    {
        amrex::ParmParse pp("CoriolisForcing");
        pp.add("latitude", 45.0);
    }

    // Needed for initial conditions
    {
        amrex::ParmParse pp("incflo");
        pp.add("use_godunov", 1);
        amrex::Vector<std::string> physics{"ABL"};
        pp.addarr("physics", physics);

        pp.add("density", 1.0); // Density
        amrex::Vector<amrex::Real> vel{{20.0, 10.0, 0.0}};
        pp.addarr("velocity", vel);
        amrex::Vector<amrex::Real> grav{{0.0, 0.0, -9.81}};
        pp.addarr("gravity", grav);
    }

    // Adjust computational domain to be more like ABL mesh in the z direction
    {
        amrex::ParmParse pp("amr");
        amrex::Vector<int> ncell{{8, 8, 64}};
        pp.addarr("n_cell", ncell);
    }
    {
        amrex::ParmParse pp("geometry");
        amrex::Vector<amrex::Real> probhi{{120.0, 120.0, 1000.0}};
        pp.addarr("prob_hi", probhi);
    }
}

} // namespace amr_wind_tests::utils
