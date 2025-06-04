#include "amr-wind/utilities/raster_ops.H"
#include "amr-wind/utilities/linear_interpolation.H"

#include <gtest/gtest.h>
#include <iostream>

using amr_wind::utils::RasterASC;

double simple_bilinear(
    double xll,
    double yll,
    double dx,
    int nx,
    int ny,
    const std::vector<double>& vals,
    double x,
    double y)
{
    double fx = (x - xll) / dx;
    double fy = (y - yll) / dx;
    int i = static_cast<int>(fx);
    int j = static_cast<int>(fy);

    // Clamp indices (for safety)
    if (i < 0) i = 0;
    if (j < 0) j = 0;
    if (i >= nx - 1) i = nx - 2;
    if (j >= ny - 1) j = ny - 2;

    double ddx = fx - i;
    double ddy = fy - j;

    int idx00 = j * nx + i;
    int idx10 = j * nx + (i + 1);
    int idx01 = (j + 1) * nx + i;
    int idx11 = (j + 1) * nx + (i + 1);

    double v00 = vals[idx00];
    double v10 = vals[idx10];
    double v01 = vals[idx01];
    double v11 = vals[idx11];

    double v0 = v00 * (1.0 - ddx) + v10 * ddx;
    double v1 = v01 * (1.0 - ddx) + v11 * ddx;
    return v0 * (1.0 - ddy) + v1 * ddy;
}

TEST(RasterASC, Interpolation)
{
    RasterASC raster;

    // Set up a 2x2 raster for testing
    amrex::Vector<amrex::Real> vals = {1.0, 2.0, 3.0, 4.0};
    raster.from_data(2, 2, 0.0, 0.0, 1.0, -9999.0, vals);

    // Test interpolation at grid points
    EXPECT_DOUBLE_EQ(raster.interp(0.0, 0.0), 1.0);
    EXPECT_DOUBLE_EQ(raster.interp(1.0, 0.0), 2.0);
    EXPECT_DOUBLE_EQ(raster.interp(0.0, 1.0), 3.0);
    EXPECT_DOUBLE_EQ(raster.interp(1.0, 1.0), 4.0);

    // Test interpolation at center (should be 2.5)
    EXPECT_NEAR(raster.interp(0.5, 0.5), 2.5, 1e-12);

    // Test nodata handling
    vals[0] = -9999.0;
    raster.from_data(2, 2, 0.0, 0.0, 1.0, -9999.0, vals);
    EXPECT_DOUBLE_EQ(raster.interp(0.0, 0.0), raster.nodata_value());
}

TEST(RasterASC, ReadAndInterp)
{
    RasterASC raster;
    raster.read("unit_tests/utilities/raster.asc");

    std::vector<double> vals(
        raster.value_ptr(), raster.value_ptr() + raster.nx() * raster.ny());

    std::vector<double> xvec(raster.nx()), yvec(raster.ny());
    for (int i = 0; i < raster.nx(); ++i)
        xvec[i] = raster.x0() + i * raster.dx();
    for (int j = 0; j < raster.ny(); ++j)
        yvec[j] = raster.y0() + j * raster.dx();

    std::vector<std::tuple<double, double, double>> sample_points = {
        {18.35246582, 46.58500101, 0.44239862},
        {35.86770315, 29.33426573, 0.94196709},
        {7.64491338, 7.64373150, 0.16562960},
        {2.84609700, 42.44263114, 0.21750405},
        {29.45463558, 34.69555631, 0.95882797},
        {1.00864022, 47.52558276, 0.20468483},
        {40.78968940, 10.40461642, 0.54802177},
        {8.90942339, 8.98682098, 0.21414673},
        {14.90786991, 25.71306515, 0.87791723},
        {21.16530591, 14.27022787, 0.82893759},
        {29.98079184, 6.83519917, 0.50818880},
        {14.31508778, 17.95173032, 0.69534286},
        {22.34742923, 38.47362211, 0.74547716},
        {9.78401533, 25.19748748, 0.53652727},
        {29.02831387, 2.27607022, 0.37606392},
        {29.76969774, 8.35568206, 0.56783972},
        {3.18752806, 46.49539133, 0.22711187},
        {47.31596962, 39.61147006, 0.65539201},
        {14.92607469, 4.78593359, 0.25920084},
        {33.52741830, 21.56747219, 1.05801096},
    };
    for (size_t idx = 0; idx < sample_points.size(); ++idx) {
        const auto& [x, y, expected] = sample_points[idx];
        double actual = raster.interp(x, y);
        double bilinear = amr_wind::interp::bilinear(xvec, yvec, vals, x, y);

        std::cout << "Point [" << idx << "] interp(" << x << ", " << y
                  << ") = " << actual << ", bilinear = " << bilinear
                  << ", expected = " << expected << std::endl;
        EXPECT_NEAR(actual, expected, 1e-8);
    }
}