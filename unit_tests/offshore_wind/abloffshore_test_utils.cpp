#include "abloffshore_test_utils.H"
#include "aw_test_utils/MeshTest.H"

namespace amr_wind_tests {

void ABLOffshoreMeshTest::populate_parameters()
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

    // Transport
    {
        amrex::ParmParse pp("transport");
        pp.add("reference_temperature", 300.0);
    }

    // ABL Forcing
    {
        amrex::ParmParse pp("ABLForcing");
        pp.add("abl_forcing_height", 90.0);
        pp.add("abl_forcing_off_height", 10.0);
        pp.add("abl_forcing_ramp_height", 30.0);
    }

    // Geostrophic Forcing
    {
        amrex::ParmParse pp("GeostrophicForcing");
        amrex::Vector<amrex::Real> gwind{{10.0, 6.0, 0.0}};
        pp.addarr("geostrophic_wind", gwind);
        pp.add("wind_forcing_off_height", 10.0);
        pp.add("wind_forcing_ramp_height", 30.0);
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
        amrex::Vector<std::string> physics{"MultiPhase", "ABL"};
        pp.addarr("physics", physics);

        pp.add("density", 1.0); // Density
        amrex::Vector<amrex::Real> vel{{20.0, 10.0, 0.0}};
        pp.addarr("velocity", vel);
        amrex::Vector<amrex::Real> grav{{0.0, 0.0, -9.81}};
        pp.addarr("gravity", grav);
    }

    // Water level
    {
        amrex::ParmParse pp("MultiPhase");
        // half cell height
        pp.add("water_level", 0.5 * 1000.0 / 128.0);
    }

    // Adjust computational domain to be more like ABL mesh in the z direction
    // with a negative portion as well for the water to occupy
    {
        amrex::ParmParse pp("amr");
        amrex::Vector<int> ncell{{8, 8, 128}};
        pp.addarr("n_cell", ncell);
    }

    {
        amrex::ParmParse pp("geometry");
        amrex::Vector<amrex::Real> probhi{{120.0, 120.0, 500.0}};
        pp.addarr("prob_hi", probhi);
        amrex::Vector<amrex::Real> problo{{0.0, 0.0, -500.0}};
        pp.addarr("prob_lo", problo);
    }
}

} // namespace amr_wind_tests
