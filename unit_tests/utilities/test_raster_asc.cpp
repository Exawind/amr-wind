#include "amr-wind/utilities/raster_ops.H"
#include "amr-wind/utilities/linear_interpolation.H"

#include <gtest/gtest.h>
#include <iostream>

using amr_wind::utils::RasterASC;

TEST(RasterASC, SetPixelInterpretationFromString_Aliases)
{
    RasterASC raster;

    raster.set_pixel_interpretation(
        amr_wind::utils::pixel_interpretation_from_string("area"));
    EXPECT_EQ(
        raster.pixel_interpretation(),
        amr_wind::utils::PixelInterpretation::Area);

    raster.set_pixel_interpretation(
        amr_wind::utils::pixel_interpretation_from_string("cell-centered"));
    EXPECT_EQ(
        raster.pixel_interpretation(),
        amr_wind::utils::PixelInterpretation::Area);

    raster.set_pixel_interpretation(
        amr_wind::utils::pixel_interpretation_from_string("cellcentered"));
    EXPECT_EQ(
        raster.pixel_interpretation(),
        amr_wind::utils::PixelInterpretation::Area);

    raster.set_pixel_interpretation(
        amr_wind::utils::pixel_interpretation_from_string("point"));
    EXPECT_EQ(
        raster.pixel_interpretation(),
        amr_wind::utils::PixelInterpretation::Point);

    raster.set_pixel_interpretation(
        amr_wind::utils::pixel_interpretation_from_string("node-centered"));
    EXPECT_EQ(
        raster.pixel_interpretation(),
        amr_wind::utils::PixelInterpretation::Point);

    raster.set_pixel_interpretation(
        amr_wind::utils::pixel_interpretation_from_string("nodecentered"));
    EXPECT_EQ(
        raster.pixel_interpretation(),
        amr_wind::utils::PixelInterpretation::Point);
}

TEST(RasterASC, InvalidPixelInterpretationStringAborts)
{
    // This test expects the code to abort when given an invalid string
    EXPECT_THROW(
        {
            auto interp =
                amr_wind::utils::pixel_interpretation_from_string("Value");
            (void)interp;
        },
        std::exception);
}

TEST(RasterASC, Interpolation)
{
    RasterASC raster;

    // Set up a 2x2 raster for testing
    amrex::Vector<amrex::Real> vals = {1.0, 2.0, 3.0, 4.0};
    raster.from_data(
        2, 2, 0.0, 0.0, 1.0, -9999.0, vals,
        amr_wind::utils::PixelInterpretation::Point);

    // Test interpolation at grid points
    EXPECT_DOUBLE_EQ(raster.interp(0.0, 0.0), 1.0);
    EXPECT_DOUBLE_EQ(raster.interp(1.0, 0.0), 2.0);
    EXPECT_DOUBLE_EQ(raster.interp(0.0, 1.0), 3.0);
    EXPECT_DOUBLE_EQ(raster.interp(1.0, 1.0), 4.0);

    // Test interpolation at center (should be 2.5)
    EXPECT_NEAR(raster.interp(0.5, 0.5), 2.5, 1e-12);

    // Test nodata handling
    vals[0] = -9999.0;
    raster.from_data(
        2, 2, 0.0, 0.0, 1.0, -9999.0, vals,
        amr_wind::utils::pixel_interpretation_from_string("point"));
    EXPECT_DOUBLE_EQ(raster.interp(0.0, 0.0), raster.nodata_value());
}

TEST(RasterASC, ReadAndInterp)
{
    // Check if the file exists
    std::ifstream infile("unit_tests/utilities/raster.asc");
    if (!infile.good()) {
        GTEST_SKIP() << "raster.asc not found, skipping test.";
    }

    RasterASC raster;
    raster.read("unit_tests/utilities/raster.asc");

    amrex::Vector<amrex::Real> vals(
        raster.value_ptr(),
        raster.value_ptr() + static_cast<std::size_t>(raster.nx()) *
                                 static_cast<std::size_t>(raster.ny()));

    amrex::Vector<amrex::Real> xvec(raster.nx()), yvec(raster.ny());
    for (int i = 0; i < raster.nx(); ++i) {
        xvec[i] = raster.x0() + i * raster.dx();
    }
    for (int j = 0; j < raster.ny(); ++j) {
        yvec[j] = raster.y0() + j * raster.dx();
    }
    const auto nx = raster.nx();
    const auto ny = raster.ny();
    amrex::Vector<amrex::Real> vals_colmajor(
        static_cast<std::size_t>(nx) * static_cast<std::size_t>(ny));
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            vals_colmajor
                [static_cast<std::size_t>(i) * static_cast<std::size_t>(ny) +
                 static_cast<std::size_t>(j)] = vals
                    [static_cast<std::size_t>(j) *
                         static_cast<std::size_t>(nx) +
                     static_cast<std::size_t>(i)];
        }
    }

    amrex::Vector<
        std::tuple<amrex::Real, amrex::Real, amrex::Real, amrex::Real>>
        sample_points = {
            {18.35246582, 56.09214408, 0.37614633, 0.37233701},
            {35.86770315, 35.32085057, 0.76834714, 0.79870701},
            {7.64491338, 9.20367670, 0.18235972, 0.16505601},
            {2.84609700, 51.10439260, 0.20160817, 0.19237600},
            {29.45463558, 41.77628209, 0.63283002, 0.65832299},
            {1.00864022, 57.22468128, 0.20113749, 0.20032200},
            {40.78968940, 12.52800753, 0.57948968, 0.57772398},
            {8.90942339, 10.82086608, 0.24168717, 0.20292400},
            {14.90786991, 30.96062947, 0.76456552, 0.73534501},
            {21.16530591, 17.18251927, 0.99193907, 0.97666001},
            {29.98079184, 8.23013778, 0.55662496, 0.54832298},
            {14.31508778, 21.61534875, 0.78565062, 0.75199401},
            {22.34742923, 46.32538172, 0.46843597, 0.46941999},
            {9.78401533, 30.33983187, 0.48460924, 0.44859299},
            {29.02831387, 2.74057435, 0.38430209, 0.36856499},
            {29.76969774, 10.06092330, 0.63800644, 0.63917899},
            {3.18752806, 55.98424670, 0.21921830, 0.21427999},
            {47.31596962, 47.69544354, 0.63721744, 0.63299602},
            {14.92607469, 5.76265473, 0.27800605, 0.24274500},
            {33.52741830, 25.96899713, 1.07977482, 1.10926795},
        };
    for (const auto& [x, y, expected, expected_nearest] : sample_points) {
        raster.set_pixel_interpretation(
            amr_wind::utils::PixelInterpretation::Point);
        amrex::Real actual = raster.interp(x, y);
        amrex::Real bilinear =
            amr_wind::interp::bilinear(xvec, yvec, vals_colmajor, x, y);

        EXPECT_NEAR(actual, expected, 1e-8);
        EXPECT_NEAR(actual, bilinear, 1e-8);

        // This is the standard for elevation, roughness, etc ...
        raster.set_pixel_interpretation(
            amr_wind::utils::PixelInterpretation::Area);
        amrex::Real actual_nearest = raster.find_nearest(x, y);
        EXPECT_NEAR(actual_nearest, expected_nearest, 1e-6);
    }
}