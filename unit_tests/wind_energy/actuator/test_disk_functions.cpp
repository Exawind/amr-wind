#include "gtest/gtest.h"
#include "amr-wind/wind_energy/actuator/disk/disk_ops.H"

namespace amr_wind::actuator::disk {

struct AreaComputerData
{
    amrex::Real radius;
    int num_points_r;
    int num_points_theta;
};

class TestAreaComputer : public ::testing::TestWithParam<AreaComputerData>
{
public:
    struct PrintParamNamesToString
    {
        template <class ParamType>
        std::string operator()(const testing::TestParamInfo<ParamType>& info)
        {
            std::stringstream name;
            const auto data = static_cast<AreaComputerData>(info.param);
            name << "R" << data.radius;
            name << "nR" << data.num_points_r;
            name << "nT" << data.num_points_theta;
            return name.str();
        }
    };

protected:
    void SetUp() final
    {
        const auto params = GetParam();
        computer = std::make_unique<ops::base::AreaComputer>(
            params.radius, params.num_points_r, params.num_points_theta);
        area = M_PI * params.radius * params.radius;
    }
    std::unique_ptr<ops::base::AreaComputer> computer;
    amrex::Real area;
};

// cppcheck-suppress uninitDerivedMemberVar
TEST_P(TestAreaComputer, area_matches)
{
    const auto params = GetParam();
    amrex::Real area_computed = 0.0;

    for (int i = 0; i < params.num_points_r; ++i) {
        for (int j = 0; j < params.num_points_theta; ++j) {
            area_computed += computer->area_section(i);
        }
    }
    EXPECT_NEAR(area, area_computed, 1.e-12 * area);
}

// cppcheck-suppress uninitDerivedMemberVar
TEST_P(TestAreaComputer, weight_sums_to_one)
{
    const auto params = GetParam();
    amrex::Real weight_computed = 0.0;

    for (int i = 0; i < params.num_points_r; ++i) {
        for (int j = 0; j < params.num_points_theta; ++j) {
            weight_computed += computer->weight(i);
        }
    }
    EXPECT_DOUBLE_EQ(1.0, weight_computed);
}

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

INSTANTIATE_TEST_SUITE_P(
    ManyAreas,
    TestAreaComputer,
    testing::Values(
        AreaComputerData{1.0, 3, 1},
        AreaComputerData{1.0, 3, 2},
        AreaComputerData{127.0, 10, 3},
        AreaComputerData{30.0, 5, 30}),
    TestAreaComputer::PrintParamNamesToString());
} // namespace amr_wind::actuator::disk
