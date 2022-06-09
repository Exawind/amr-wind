#include "gtest/gtest.h"
#include "amr-wind/wind_energy/actuator/disk/disk_ops.H"

namespace amr_wind {
namespace actuator {
namespace disk {

class TestComputeDiskPoints : public ::testing::TestWithParam<vs::Vector>
{};

TEST_P(TestComputeDiskPoints, disk_translation_x_aligned)
{
    DiskBaseData data;
    data.center = GetParam();
    data.normal_vec = {1.0, 0.0, 0.0};
    data.num_vel_pts_r = 2;
    data.num_vel_pts_t = 2;
    data.diameter = 2.0;

    const amrex::Real dr = data.diameter * 0.5 / data.num_vel_pts_r;

    VecList points(4);
    VecList gold_points(4);
    const amrex::Real eps = 3.0 * vs::DTraits<amrex::Real>::eps();
    gold_points[0] = vs::Vector{0.0, 0.5 * dr, 0.0} + data.center;
    gold_points[1] = vs::Vector{0.0, -0.5 * dr, 0.0} + data.center;
    gold_points[2] = vs::Vector{0.0, 1.5 * dr, 0.0} + data.center;
    gold_points[3] = vs::Vector{0.0, -1.5 * dr, 0.0} + data.center;
    ops::base::compute_disk_points(data, points, data.normal_vec, 0, 0);
    for (int i = 0; i < points.size(); ++i) {
        EXPECT_NEAR(gold_points[i].x(), points[i].x(), eps);
        EXPECT_NEAR(gold_points[i].y(), points[i].y(), eps);
        EXPECT_NEAR(gold_points[i].z(), points[i].z(), eps);
    }
}

INSTANTIATE_TEST_SUITE_P(
    ChangeCenters,
    TestComputeDiskPoints,
    testing::Values(
        vs::Vector{0, 0, 0},
        vs::Vector{1, 0, 0},
        vs::Vector{0, 1, 0},
        vs::Vector{0, 0, 1},
        vs::Vector{1, 1, 1},
        vs::Vector{-1, 0, 0},
        vs::Vector{-1, -1, -1},
        vs::Vector{1000.0, 30.0, 5.0}));
} // namespace disk
} // namespace actuator
} // namespace amr_wind