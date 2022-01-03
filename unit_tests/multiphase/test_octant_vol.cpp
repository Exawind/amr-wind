#include "aw_test_utils/MeshTest.H"
#include "aw_test_utils/iter_tools.H"
#include "aw_test_utils/test_utils.H"
#include "amr-wind/equation_systems/vof/volume_fractions.H"

namespace amr_wind_tests {

class VOFInterpTest : public AmrexTest {};

TEST_F(VOFInterpTest, Octants1D)
{
    const amrex::Real tol = 1e-8;

    // Cell has VOF = 0.5, interface plane cuts vertically with liquid on left
    // side, only varies in x

    amrex::Real m1 = 1.0;
    amrex::Real m2 = 0.0;
    amrex::Real m3 = 0.0;
    amrex::Real alpha = 0.5;

    auto vol_oct = amr_wind::multiphase::cut_volume_oct(m1,m2,m3,alpha);

    // Left, bottom, back
    EXPECT_NEAR(vol_oct[0], 1.0, tol);
    // Right, bottom, back
    EXPECT_NEAR(vol_oct[1], 0.0, tol);
    // Left, top, back
    EXPECT_NEAR(vol_oct[2], 1.0, tol);
    // Right, top, back
    EXPECT_NEAR(vol_oct[3], 0.0, tol);
    // Left, bottom, front
    EXPECT_NEAR(vol_oct[4], 1.0, tol);
    // Right, bottom, front
    EXPECT_NEAR(vol_oct[5], 0.0, tol);
    // Left, top, front
    EXPECT_NEAR(vol_oct[6], 1.0, tol);
    // Right, top, front
    EXPECT_NEAR(vol_oct[7], 0.0, tol);

}

TEST_F(VOFInterpTest, Octants2D)
{
    const amrex::Real tol = 1e-8;

    // Cell has VOF = 0.5, interface plane cuts diagonally through cell with
    // liquid on the left side, no variation in z

    amrex::Real m1 = 0.5;
    amrex::Real m2 = 0.5;
    amrex::Real m3 = 0.0;
    amrex::Real alpha = 0.5;

    auto vol_oct = amr_wind::multiphase::cut_volume_oct(m1,m2,m3,alpha);

    // Left, bottom, back
    EXPECT_NEAR(vol_oct[0], 1.0, tol);
    // Right, bottom, back
    EXPECT_NEAR(vol_oct[1], 0.5, tol);
    // Left, top, back
    EXPECT_NEAR(vol_oct[2], 0.5, tol);
    // Right, top, back
    EXPECT_NEAR(vol_oct[3], 0.0, tol);
    // Left, bottom, front
    EXPECT_NEAR(vol_oct[4], 1.0, tol);
    // Right, bottom, front
    EXPECT_NEAR(vol_oct[5], 0.5, tol);
    // Left, top, front
    EXPECT_NEAR(vol_oct[6], 0.5, tol);
    // Right, top, front
    EXPECT_NEAR(vol_oct[7], 0.0, tol);

}

TEST_F(VOFInterpTest, Octants3D)
{
    const amrex::Real tol = 1e-8;

    // Cell has VOF = 0.5, interface plane cuts diagonally through cell with
    // liquid on the left side, intersects with (0, 1, 0) and (1, 0, 1)

    amrex::Real m1 = 1.0 / 4.0;
    amrex::Real m2 = 2.0 / 4.0;
    amrex::Real m3 = 1.0 / 4.0;
    amrex::Real alpha = 0.5;

    auto vol_oct = amr_wind::multiphase::cut_volume_oct(m1,m2,m3,alpha);

    // Left, bottom, back
    EXPECT_NEAR(vol_oct[0], 1.0, tol);
    // Right, bottom, back
    //EXPECT_NEAR(vol_oct[1], 0.5, tol);
    // Left, top, back
    EXPECT_NEAR(vol_oct[2], 0.5, tol);
    // Right, top, back
    //EXPECT_NEAR(vol_oct[3], 0.0, tol);
    // Left, bottom, front
    //EXPECT_NEAR(vol_oct[4], 1.0, tol);
    // Right, bottom, front
    EXPECT_NEAR(vol_oct[5], 0.5, tol);
    // Left, top, front
    //EXPECT_NEAR(vol_oct[6], 0.5, tol);
    // Right, top, front
    EXPECT_NEAR(vol_oct[7], 0.0, tol);

}

} // namespace amr_wind_tests
