#include "aw_test_utils/MeshTest.H"
#include "aw_test_utils/iter_tools.H"
#include "aw_test_utils/test_utils.H"
#include "amr-wind/ocean_waves/utils/wave_utils_K.H"

namespace amr_wind_tests {

class WaveUtilsTest : public MeshTest
{};

TEST_F(WaveUtilsTest, free_surface_to_vof)
{
    constexpr amrex::Real tol = 1e-12;
    constexpr amrex::Real zsl = 1.0;
    constexpr amrex::Real dz = 0.5;

    const amrex::Real vof_fully_above =
        amr_wind::ocean_waves::utils::free_surface_to_vof(zsl, 100., dz);
    const amrex::Real vof_fully_below =
        amr_wind::ocean_waves::utils::free_surface_to_vof(zsl, -100., dz);
    const amrex::Real vof_at_level =
        amr_wind::ocean_waves::utils::free_surface_to_vof(zsl, zsl, dz);
    const amrex::Real vof_half_above =
        amr_wind::ocean_waves::utils::free_surface_to_vof(
            zsl, zsl + 0.25 * dz, dz);
    const amrex::Real vof_half_below =
        amr_wind::ocean_waves::utils::free_surface_to_vof(
            zsl, zsl - 0.25 * dz, dz);
    const amrex::Real vof_exactly_above =
        amr_wind::ocean_waves::utils::free_surface_to_vof(
            zsl, zsl + 0.5 * dz, dz);
    const amrex::Real vof_exactly_below =
        amr_wind::ocean_waves::utils::free_surface_to_vof(
            zsl, zsl - 0.5 * dz, dz);

    EXPECT_NEAR(vof_fully_above, 0.0, tol);
    EXPECT_NEAR(vof_fully_below, 1.0, tol);
    EXPECT_NEAR(vof_at_level, 0.5, tol);
    EXPECT_NEAR(vof_half_above, 0.25, tol);
    EXPECT_NEAR(vof_half_below, 0.75, tol);
    EXPECT_NEAR(vof_exactly_above, 0.0, tol);
    EXPECT_NEAR(vof_exactly_below, 1.0, tol);
}

TEST_F(WaveUtilsTest, gamma_generate)
{
    constexpr amrex::Real tol = 1e-12;
    constexpr amrex::Real zone_length = 1.5;

    const amrex::Real gamma_past_left =
        amr_wind::ocean_waves::utils::gamma_generate(-1.0, zone_length);
    const amrex::Real gamma_left =
        amr_wind::ocean_waves::utils::gamma_generate(0.0, zone_length);
    const amrex::Real gamma_middle =
        amr_wind::ocean_waves::utils::gamma_generate(
            0.5 * zone_length, zone_length);
    const amrex::Real gamma_right =
        amr_wind::ocean_waves::utils::gamma_generate(zone_length, zone_length);
    const amrex::Real gamma_past_right =
        amr_wind::ocean_waves::utils::gamma_generate(
            zone_length + 1.0, zone_length);

    const amrex::Real gamma_middle_gold =
        1.0 - (std::exp(std::pow(0.5, 3.5)) - 1.0) / (std::exp(1.0) - 1.0);

    EXPECT_NEAR(gamma_past_left, 0.0, tol);
    EXPECT_NEAR(gamma_left, 0.0, tol);
    EXPECT_NEAR(gamma_middle, gamma_middle_gold, tol);
    EXPECT_NEAR(gamma_right, 1.0, tol);
    EXPECT_NEAR(gamma_past_right, 1.0, tol);
}

TEST_F(WaveUtilsTest, gamma_absorb)
{
    constexpr amrex::Real tol = 1e-12;
    constexpr amrex::Real zone_length = 1.5;

    const amrex::Real gamma_past_left =
        amr_wind::ocean_waves::utils::gamma_absorb(-1.0, zone_length, 1.0);
    const amrex::Real gamma_left =
        amr_wind::ocean_waves::utils::gamma_absorb(0.0, zone_length, 1.0);
    const amrex::Real gamma_middle = amr_wind::ocean_waves::utils::gamma_absorb(
        0.5 * zone_length, zone_length, 1.0);
    const amrex::Real gamma_right = amr_wind::ocean_waves::utils::gamma_absorb(
        zone_length, zone_length, 1.0);
    const amrex::Real gamma_past_right =
        amr_wind::ocean_waves::utils::gamma_absorb(
            zone_length + 1.0, zone_length, 1.0);

    const amrex::Real gamma_right_length_factor2 =
        amr_wind::ocean_waves::utils::gamma_absorb(
            zone_length, zone_length, 2.0);

    const amrex::Real gamma_middle_gold =
        1.0 - (std::exp(std::pow(0.5, 3.5)) - 1.0) / (std::exp(1.0) - 1.0);

    EXPECT_NEAR(gamma_past_left, 1.0, tol);
    EXPECT_NEAR(gamma_left, 1.0, tol);
    EXPECT_NEAR(gamma_middle, gamma_middle_gold, tol);
    EXPECT_NEAR(gamma_right_length_factor2, gamma_middle_gold, tol);
    EXPECT_NEAR(gamma_right, 0.0, tol);
    EXPECT_NEAR(gamma_past_right, 0.0, tol);
}

TEST_F(WaveUtilsTest, ramp)
{
    constexpr amrex::Real tol = 1e-12;
    constexpr amrex::Real near_tol = 1e-5;
    constexpr amrex::Real ramp_period = 1.5;

    const amrex::Real f_ramp_begin =
        amr_wind::ocean_waves::utils::ramp(0.0, ramp_period);
    const amrex::Real f_ramp_middle =
        amr_wind::ocean_waves::utils::ramp(0.5 * ramp_period, ramp_period);
    const amrex::Real f_ramp_near_end =
        amr_wind::ocean_waves::utils::ramp(ramp_period - near_tol, ramp_period);
    const amrex::Real f_ramp_end =
        amr_wind::ocean_waves::utils::ramp(ramp_period, ramp_period);
    const amrex::Real f_ramp_past =
        amr_wind::ocean_waves::utils::ramp(ramp_period + 1.0, ramp_period);

    const amrex::Real f_ramp_middle_gold = 0.5 - std::sin(M_PI * 0.5) / M_PI;

    EXPECT_NEAR(f_ramp_begin, 0.0, tol);
    EXPECT_NEAR(f_ramp_middle, f_ramp_middle_gold, tol);
    EXPECT_NEAR(f_ramp_near_end, 1.0, 2.0 * near_tol);
    EXPECT_NEAR(f_ramp_end, 1.0, tol);
    EXPECT_NEAR(f_ramp_past, 1.0, tol);
}

TEST_F(WaveUtilsTest, combine_linear)
{
    constexpr amrex::Real tol = 1e-12;
    constexpr amrex::Real Gamma = 0.7;
    constexpr amrex::Real target = 2.0;
    constexpr amrex::Real current = 1.5;

    const amrex::Real result =
        amr_wind::ocean_waves::utils::combine_linear(Gamma, target, current);

    EXPECT_NEAR(result, (1.0 - Gamma) * target + Gamma * current, tol);
}

TEST_F(WaveUtilsTest, harmonize_profiles)
{
    constexpr amrex::Real tol = 1e-12;
    constexpr amrex::Real problo_x = -1.0;
    constexpr amrex::Real probhi_x = 1.0;
    constexpr amrex::Real gen_length = 0.2;
    constexpr amrex::Real beach_length = 0.4;

    const amr_wind::ocean_waves::utils::WaveVec left{1.0, 2.0, 3.0, 0.1};
    const amr_wind::ocean_waves::utils::WaveVec bulk{0.0, 1.1, 2.0, 0.0};
    const amr_wind::ocean_waves::utils::WaveVec right{-1.0, 1.0, 0.5, -0.1};

    amrex::Real x = -0.9;
    auto result = amr_wind::ocean_waves::utils::harmonize_profiles_1d(
        x, problo_x, gen_length, probhi_x, beach_length, left, bulk, right);
    for (int n = 0; n < 4; ++n) {
        EXPECT_NEAR(left[n], result[n], tol);
    }

    x = -0.75;
    result = amr_wind::ocean_waves::utils::harmonize_profiles_1d(
        x, problo_x, gen_length, probhi_x, beach_length, left, bulk, right);
    const amrex::Real Gamma_l = amr_wind::ocean_waves::utils::gamma_generate(
        x - (problo_x + gen_length), 0.5 * gen_length);
    for (int n = 0; n < 4; ++n) {
        EXPECT_NEAR(
            (1.0 - Gamma_l) * left[n] + Gamma_l * bulk[n], result[n], tol);
    }

    x = 0.;
    result = amr_wind::ocean_waves::utils::harmonize_profiles_1d(
        x, problo_x, gen_length, probhi_x, beach_length, left, bulk, right);
    for (int n = 0; n < 4; ++n) {
        EXPECT_NEAR(bulk[n], result[n], tol);
    }

    x = 0.5;
    result = amr_wind::ocean_waves::utils::harmonize_profiles_1d(
        x, problo_x, gen_length, probhi_x, beach_length, left, bulk, right);
    const amrex::Real Gamma_r = amr_wind::ocean_waves::utils::gamma_absorb(
        x - (probhi_x - beach_length) + 0.5 * beach_length, 0.5 * beach_length,
        1.0);
    for (int n = 0; n < 4; ++n) {
        EXPECT_NEAR(
            (1.0 - Gamma_r) * right[n] + Gamma_r * bulk[n], result[n], tol);
    }

    x = 0.8;
    result = amr_wind::ocean_waves::utils::harmonize_profiles_1d(
        x, problo_x, gen_length, probhi_x, beach_length, left, bulk, right);
    for (int n = 0; n < 4; ++n) {
        EXPECT_NEAR(right[n], result[n], tol);
    }
}

} // namespace amr_wind_tests