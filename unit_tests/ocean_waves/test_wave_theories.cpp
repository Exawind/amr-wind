#include "aw_test_utils/MeshTest.H"
#include "aw_test_utils/iter_tools.H"
#include "aw_test_utils/test_utils.H"
#include "amr-wind/ocean_waves/relaxation_zones/stokes_waves_K.H"

namespace amr_wind_tests {

class WaveTheoriesTest : public MeshTest
{};

TEST_F(WaveTheoriesTest, StokesWaves)
{
    amrex::Real CoeffTol = 1e-4;

    amrex::Real wavenumber = 2.0;
    amrex::Real waterdepth = 0.376991;
    int StokesOrder = 5;
    amrex::Real c0, a11, a22, b22, c2, d2, e2, a31, a33, b31, a42;
    amrex::Real a44, b42, b44, c4, d4, e4, a51, a53, a55, b53, b55;

    amr_wind::ocean_waves::relaxation_zones::stokes_coefficients(
        StokesOrder, wavenumber, waterdepth, c0, a11, a22, b22, c2, d2, e2, a31,
        a33, b31, a42, a44, b42, b44, c4, d4, e4, a51, a53, a55, b53, b55);

    // Gold coefficient values taken from table 2 of
    // Fenton, J. Fifth Order Stokes Theory for Steady Waves
    // Journal of Waterway, Port, Coastal and Ocean Engineering, 1985, 111,
    // 216-234

    const amrex::Real gold_A11 = 1.208490;
    const amrex::Real gold_A22 = 0.799840;
    const amrex::Real gold_A31 = -9.105340;
    const amrex::Real gold_A33 = 0.368275;
    const amrex::Real gold_A42 = -12.196150;
    const amrex::Real gold_A44 = 0.058723;
    const amrex::Real gold_A51 = 108.46831725;
    const amrex::Real gold_A53 = -6.941756;
    const amrex::Real gold_A55 = -0.074979;
    const amrex::Real gold_B22 = 2.502414;
    const amrex::Real gold_B31 = -5.731666;
    const amrex::Real gold_B42 = -32.407508;
    const amrex::Real gold_B44 = 14.033758;
    const amrex::Real gold_B53 = -103.44536875;
    const amrex::Real gold_B55 = 37.200027;
    const amrex::Real gold_C0 = 0.798448;
    const amrex::Real gold_C2 = 1.940215;
    const amrex::Real gold_C4 = -12.970403;
    const amrex::Real gold_D2 = -0.626215;
    const amrex::Real gold_D4 = 3.257104;
    const amrex::Real gold_E2 = 1.781926;
    const amrex::Real gold_E4 = -11.573657;

    EXPECT_NEAR(gold_A11, a11, CoeffTol);
    EXPECT_NEAR(gold_A22, a22, CoeffTol);
    EXPECT_NEAR(gold_A31, a31, CoeffTol);
    EXPECT_NEAR(gold_A33, a33, CoeffTol);
    EXPECT_NEAR(gold_A42, a42, CoeffTol);
    EXPECT_NEAR(gold_A44, a44, CoeffTol);
    EXPECT_NEAR(gold_A51, a51, CoeffTol);
    EXPECT_NEAR(gold_A53, a53, CoeffTol);
    EXPECT_NEAR(gold_A55, a55, CoeffTol);
    EXPECT_NEAR(gold_B22, b22, CoeffTol);
    EXPECT_NEAR(gold_B31, b31, CoeffTol);
    EXPECT_NEAR(gold_B42, b42, CoeffTol);
    EXPECT_NEAR(gold_B44, b44, CoeffTol);
    EXPECT_NEAR(gold_B53, b53, CoeffTol);
    EXPECT_NEAR(gold_B55, b55, CoeffTol);
    EXPECT_NEAR(gold_C0, c0, CoeffTol);
    EXPECT_NEAR(gold_C2, c2, CoeffTol);
    EXPECT_NEAR(gold_C4, c4, CoeffTol);
    EXPECT_NEAR(gold_D2, d2, CoeffTol);
    EXPECT_NEAR(gold_D4, d4, CoeffTol);
    EXPECT_NEAR(gold_E2, e2, CoeffTol);
    EXPECT_NEAR(gold_E4, e4, CoeffTol);
}

TEST_F(WaveTheoriesTest, StokesWaveLength)
{

    amrex::Real wave_height = 0.2;
    amrex::Real wave_period = 2.2;
    amrex::Real water_depth = 5.0;
    int wave_order = 5;
    constexpr amrex::Real g = 9.81;
    amrex::Real lambda =
        amr_wind::ocean_waves::relaxation_zones::stokes_wave_length(
            wave_period, water_depth, wave_height, wave_order, g);

    // Relation to check is from course notes:
    // https://www.caee.utexas.edu/prof/kinnas/ce358/oenotes/kinnas_stokes11.pdf
    amrex::Real k = 2.0 * M_PI / lambda;
    const amrex::Real RHS1 = 2.0 * M_PI / (wave_period * std::sqrt(g * k));

    amrex::Real S = 1.0 / std::cosh(2.0 * k * water_depth);
    amrex::Real C = 1.0 - S;
    amrex::Real eps = k * wave_height / 2.0;

    amrex::Real C0 = std::sqrt(std::tanh(k * water_depth));
    amrex::Real C2 = C0 * (2.0 + 7.0 * S * S) / (4.0 * C * C);
    const amrex::Real C4 =
        C0 *
        (4.0 + 32.0 * S - 116.0 * std::pow(S, 2) - 400.0 * std::pow(S, 3) -
         71.0 * std::pow(S, 4) + 146.0 * std::pow(S, 5)) /
        (32.0 * std::pow(C, 5));
    const amrex::Real LHS1 = C0 + std::pow(eps, 2) * C2 + std::pow(eps, 4) * C4;
    EXPECT_NEAR(LHS1, RHS1, 1e-8);

    // Reevaluate with a new set of conditions
    wave_height = 0.05;
    wave_period = 1.5;
    water_depth = 0.9;
    wave_order = 3;
    lambda = amr_wind::ocean_waves::relaxation_zones::stokes_wave_length(
        wave_period, water_depth, wave_height, wave_order, g);

    k = 2.0 * M_PI / lambda;
    const amrex::Real RHS2 = 2.0 * M_PI / (wave_period * std::sqrt(g * k));

    S = 1.0 / std::cosh(2.0 * k * water_depth);
    C = 1.0 - S;
    eps = k * wave_height / 2.0;

    C0 = std::sqrt(std::tanh(k * water_depth));
    C2 = C0 * (2.0 + 7.0 * S * S) / (4.0 * C * C);
    const amrex::Real LHS2 = C0 + std::pow(eps, 2) * C2;
    EXPECT_NEAR(LHS2, RHS2, 1e-8);
}

} // namespace amr_wind_tests