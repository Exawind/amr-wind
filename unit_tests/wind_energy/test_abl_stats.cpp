#include "abl_test_utils.H"
#include "aw_test_utils/iter_tools.H"
#include "aw_test_utils/test_utils.H"
#include "amr-wind/equation_systems/tke/TKE.H"
#include "amr-wind/wind_energy/ABLStats.H"

namespace amr_wind_tests {

TEST_F(ABLMeshTest, stats_tke_diffusion)
{
    // This test checks the implementation of the calc_diffusion routine to see
    // if it does what it is intended to do
    constexpr double tol = 1.0e-12;
    constexpr amrex::Real val_shear = 0.5;
    constexpr amrex::Real val_buoy = 0.4;
    constexpr amrex::Real val_dissip = 0.3;
    constexpr amrex::Real val_conv = 0.2;
    constexpr amrex::Real val_tke = 1.6;
    constexpr amrex::Real val_tkeold = 1.0;
    constexpr amrex::Real dt = 0.1;
    constexpr amrex::Real expected_diff = (val_tke - val_tkeold) / dt -
                                          val_conv - val_shear - val_buoy +
                                          val_dissip;

    populate_parameters();
    initialize_mesh();

    // Register fields and eqs for the sake of a functional test
    // sim().repo().declare_field("density", 1, 3);
    sim().repo().declare_field("temperature", 1, 3);
    sim().pde_manager().register_icns();

    // Register fields used in turbulence model OneEqKsgs
    auto& shear = sim().repo().declare_field("shear_prod", 1, 1);
    auto& buoy = sim().repo().declare_field("buoy_prod", 1, 1);
    auto& dissip = sim().repo().declare_field("dissipation", 1, 1);
    auto diff = sim().repo().create_scratch_field("diffusion", 1, 1);
    // Register transport pde for tke
    auto& tke_eqn = sim().pde_manager().register_transport_pde(
        amr_wind::pde::TKE::pde_name());
    auto& tke = tke_eqn.fields().field;

    // Populate values of fields
    shear.setVal(val_shear);
    buoy.setVal(val_buoy);
    dissip.setVal(val_dissip);
    tke_eqn.fields().conv_term.setVal(val_conv);
    tke.setVal(val_tke);
    tke.state(amr_wind::FieldState::Old).setVal(val_tkeold);

    // Initialize ABL Stats
    amr_wind::ABLWallFunction wall_func(sim());
    amr_wind::ABLStats stats(sim(), wall_func, 2);

    // Calculate diffusion term
    stats.calc_tke_diffusion(*diff, buoy, shear, dissip, dt);

    // Check answer
    const auto min_val = utils::field_min(*diff, 0);
    const auto max_val = utils::field_max(*diff, 0);
    EXPECT_NEAR(expected_diff, min_val, tol);
    EXPECT_NEAR(expected_diff, max_val, tol);
}

TEST_F(ABLMeshTest, stats_energy_budget)
{
    // This test checks the assumptions behind the calc_diffusion routine

    // Set up turbulence model

    populate_parameters();
    initialize_mesh();

    // Initialize turbulence model

    // Register transport pde for tke

    // Step forward in time

    // Calculate diffusion term

    // Check assumptions in diffusion term
}

} // namespace amr_wind_tests
