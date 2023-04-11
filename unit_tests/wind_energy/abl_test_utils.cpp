#include "abl_test_utils.H"
#include "aw_test_utils/MeshTest.H"

namespace amr_wind_tests {

void ABLMeshTest::populate_parameters()
{
    // Must populate with default parameters first to signify parameters do not
    // need to be reread within initialize_mesh
    MeshTest::populate_parameters();
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

    // Initial conditions (Linear profile)
    {
        amrex::ParmParse pp("ABL");
        bool linear_profile = false;
        pp.add("linear_profile", static_cast<int>(linear_profile));

        amrex::Vector<amrex::Real> top_velocity{{20., 0.0, 0.0}};
        amrex::Vector<amrex::Real> bottom_velocity{{4.0, 0.0, 0.0}};
        pp.addarr("top_velocity", top_velocity);
        pp.addarr("bottom_velocity", bottom_velocity);
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
    }

    // Hurricane Forcing
    {
        amrex::ParmParse pp("HurricaneForcing");
        amrex::Real gradient_wind{40.0};
        amrex::Real radial_distance{40000.0};
        amrex::Real gradient_wind_radial_decay{-0.008};
        amrex::Real gradient_wind_zero_height{18000.};

        pp.add("gradient_wind", gradient_wind);
        pp.add("eyewall_radial_distance", radial_distance);
        pp.add("gradient_wind_radial_decay", gradient_wind_radial_decay);
        pp.add("gradient_wind_zero_height", gradient_wind_zero_height);
    }

    // Rayleigh damping
    {

        amrex::ParmParse pp("RayleighDamping");
        amrex::Real time_scale{40.0};
        amrex::Real damping_length{1000};
        amrex::Real full_damping_length{500};
        amrex::Vector<amrex::Real> reference_velocity{{15., 0., 0.}};

        pp.add("time_scale", time_scale);
        pp.add("damping_length", damping_length);
        pp.add("full_damping_length", full_damping_length);
        pp.addarr("reference_velocity", reference_velocity);
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
