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

} // namespace amr_wind_tests