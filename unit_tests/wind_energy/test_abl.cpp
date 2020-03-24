#include "aw_test_utils/MeshTest.H"
#include "aw_test_utils/iter_tools.H"
#include "aw_test_utils/test_utils.H"
#include "ABLFieldInit.H"

namespace amr_wind_tests {
namespace {

// Introduce namespace for the text fixture
class ABLTest : public MeshTest
{};

void populate_abl_params()
{
    amrex::ParmParse pp("abl");

    // Initial conditions (Temperature)
    amrex::Vector<amrex::Real> theights{{0.0, 650.0, 750.0, 1000.0}};
    amrex::Vector<amrex::Real> tvalues{{300.0, 300.0, 308.0, 308.75}};
    pp.add("ntemperature", theights.size());
    pp.addarr("temperature_heights", theights);
    pp.addarr("temperature_values", tvalues);
    pp.add("perturb_ref_height", 50.0);

    // Boussinesq Buoyancy
    const amrex::Real tecoeff = 1.0 / 300.0;
    pp.add("reference_temperature", 300.0);
    pp.add("thermal_expansion_coeff", tecoeff);

    // ABL Forcing
    pp.add("abl_forcing_height", 90.0);

    // Coriolis term
    pp.add("latitude", 45.0);

    pp.add("kappa", 0.41);
    pp.add("surface_roughness_z0", 0.1);

    // Needed for initial conditions
    {
        amrex::ParmParse pp("incflo");
        pp.add("ro_0", 1.0);   // Density
        pp.add("ic_u", 20.0);  // Ux
        pp.add("ic_v", 10.0);   // Uy
        pp.add("ic_w", 0.0);   // Uw
    }

    // Adjust computational domain to be more like ABL mesh in the z direction
    {
        amrex::ParmParse pp("amr");
        amrex::Vector<int> ncell{{8, 8, 64}};
        pp.addarr("n_cell", ncell);
    }
    {
        amrex::ParmParse pp("geometry");
        amrex::Vector<amrex::Real> probhi {{120.0, 120.0, 1000.0}};

        pp.addarr("prob_hi", probhi);
    }
}

} // namespace

TEST_F(ABLTest, abl_initialization)
{
    populate_parameters();
    populate_abl_params();

    // For the first test disable velocity perturbations
    {
        amrex::ParmParse pp("abl");
        bool perturb_velocity = false;
        pp.add("perturb_velocity", static_cast<int>(perturb_velocity));
    }

    initialize_mesh();
    auto velocity = mesh().declare_field("velocity", 3, 0);
    auto density = mesh().declare_field("density");
    auto temperature = mesh().declare_field("temperature");

    amr_wind::ABLFieldInit ablinit;
    run_algorithm(
        mesh().num_levels(), density,
        [&](const int lev, const amrex::MFIter& mfi) {
            auto vel = velocity[lev]->array(mfi);
            auto rho = density[lev]->array(mfi);
            auto theta = temperature[lev]->array(mfi);

            const auto& bx = mfi.validbox();
            ablinit(bx, mesh().Geom(lev), vel, rho, theta);
        });

    const int nlevels = mesh().num_levels();
    const amrex::Real tol = 1.0e-12;

    // Test temperature
    {
        const amrex::Real dz = mesh().Geom(0).CellSize(2);
        const amrex::Real min_temp_gold = 300.0;
        const amrex::Real max_temp_gold = 308.75 - (308.75 - 308.0) / 250.0 * (0.5 * dz);
        amrex::Real min_temp, max_temp;
        utils::field_minmax(nlevels, temperature, min_temp, max_temp);
        EXPECT_NEAR(min_temp, min_temp_gold, tol);
        EXPECT_NEAR(max_temp, max_temp_gold, tol);
    }

    // Test density
    {
        amrex::Real min_rho, max_rho;
        utils::field_minmax(nlevels, density, min_rho, max_rho);
        EXPECT_NEAR(min_rho, 1.0, tol);
        EXPECT_NEAR(max_rho, 1.0, tol);
    }

    // Test velocity
    {
        amrex::Vector<amrex::Real> min_vel(3), max_vel(3);
        utils::field_minmax(nlevels, velocity, min_vel, max_vel);
        EXPECT_NEAR(min_vel[0], 20.0, tol);
        EXPECT_NEAR(min_vel[1], 10.0, tol);
        EXPECT_NEAR(min_vel[2],  0.0, tol);
        EXPECT_NEAR(max_vel[0], 20.0, tol);
        EXPECT_NEAR(max_vel[1], 10.0, tol);
        EXPECT_NEAR(max_vel[2],  0.0, tol);
    }
}

} // namespace amr_wind_tests
