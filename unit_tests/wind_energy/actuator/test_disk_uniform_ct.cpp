#include "aw_test_utils/AmrexTest.H"
#include "test_act_utils.H"

#include "amr-wind/wind_energy/actuator/disk/uniform_ct_ops.H"
#include "AMReX_Exception.H"

namespace amr_wind_tests {

namespace {
class UniformCtTest : public AmrexTest
{
protected:
    void SetUp() override
    {
        AmrexTest::SetUp();
        {
            amrex::ParmParse pp("Actuator");
            amrex::ParmParse pu("Actuator.UniformCtDisk");
            pu.add("num_force_points", 3);
            pu.add("epsilon", 1);
            pu.add("rotor_diameter", 1);
            std::vector<double> ct{1};
            pu.addarr("thrust_coeff", ct);
        }
        {
            amrex::ParmParse pp("Coriolis");
        }
    }
};
} // namespace

namespace act = amr_wind::actuator;
namespace ops = act::ops;
namespace vs = amr_wind::vs;

TEST_F(UniformCtTest, compute_vecs_from_yaw)
{
    act::UniformCtData meta;
    amrex::ParmParse pp("Actuator.UniformCtDisk");
    pp.add("yaw", 90.0);
    pp.add("sample_yaw", 45.0);
    act::utils::ActParser ap("Actuator.UniformCtDisk", "Actuator");
    ASSERT_TRUE(ap.contains("yaw"));
    ops::optional_parameters(meta, ap);
    {
        const vs::Vector gold_vec = {1, 0, 0};
        for (int i = 0; i < 3; i++) {
            EXPECT_NEAR(gold_vec[i], meta.normal_vec[i], 1e-12) << i;
        }
    }
    {
        vs::Vector gold_vec = {1, 1, 0};
        gold_vec.normalize();

        for (int i = 0; i < 3; i++) {
            EXPECT_NEAR(gold_vec[i], meta.sample_vec[i], 1e-12) << i;
        }
    }
}

TEST_F(UniformCtTest, compute_vecs_from_tilt)
{
    act::UniformCtData meta;
    amrex::ParmParse pp("Actuator.UniformCtDisk");
    pp.add("yaw", 270.0);
    pp.add("tilt", -90.0);
    pp.add("sample_yaw", -315.0);
    pp.add("sample_tilt", -45.0);
    act::utils::ActParser ap("Actuator.UniformCtDisk", "Actuator");
    ASSERT_TRUE(ap.contains("tilt"));
    ops::optional_parameters(meta, ap);
    {
        const vs::Vector gold_vec = {0, 0, 1};
        for (int i = 0; i < 3; i++) {
            EXPECT_NEAR(gold_vec[i], meta.normal_vec[i], 1e-12) << i;
        }
    }
    {
        vs::Vector gold_vec = {0.5, 0.5, std::sqrt(0.5)};
        gold_vec.normalize();

        for (int i = 0; i < 3; i++) {
            EXPECT_NEAR(gold_vec[i], meta.sample_vec[i], 1e-12) << i;
        }
    }
    // check that we throw because we've created a bad normal with this much
    // tilt
    EXPECT_THROW(
        ops::compute_and_normalize_coplanar_vector(meta), amrex::RuntimeError);
}

TEST_F(UniformCtTest, compute_vecs_with_different_north)
{
    act::UniformCtData meta;
    const std::vector<double> north{1, 0, 0};
    const std::vector<double> east{0, -1, 0};
    {
        amrex::ParmParse pp("Coriolis.Forcing");
        pp.addarr("north_vector", north);
        pp.addarr("east_vector", east);
    }
    {
        amrex::ParmParse pp("Actuator.UniformCtDisk");
    }
    act::utils::ActParser ap("Actuator.UniformCtDisk", "Actuator");
    act::utils::ActParser cp("Coriolis.Forcing", "Coriolis");
    ASSERT_TRUE(cp.contains("north_vector"));
    ops::optional_parameters(meta, ap);
    {
        const auto& gold_vec = north;
        for (int i = 0; i < 3; i++) {
            EXPECT_NEAR(gold_vec[i], meta.normal_vec[i], 1e-12) << i;
            EXPECT_NEAR(gold_vec[i], meta.sample_vec[i], 1e-12) << i;
        }
    }
}

TEST_F(UniformCtTest, compute_vecs_from_yaw_and_tilt_with_different_north)
{
    act::UniformCtData meta;
    const std::vector<double> north{1, 0, 0};
    const std::vector<double> east{0, -1, 0};
    {
        amrex::ParmParse pp("Coriolis.Forcing");
        pp.addarr("north_vector", north);
        pp.addarr("east_vector", east);
    }
    {
        amrex::ParmParse pp("Actuator.UniformCtDisk");
        pp.add("yaw", 90.0);
        pp.add("sample_yaw", 45.0);
        pp.add("sample_tilt", -45.0);
    }
    act::utils::ActParser ap("Actuator.UniformCtDisk", "Actuator");
    ops::optional_parameters(meta, ap);
    {
        const auto& gold_vec = east;
        for (int i = 0; i < 3; i++) {
            EXPECT_NEAR(gold_vec[i], meta.normal_vec[i], 1e-12) << i;
        }
    }
    {
        vs::Vector gold_vec = {0.5, -0.5, std::sqrt(0.5)};
        gold_vec.normalize();

        for (int i = 0; i < 3; i++) {
            EXPECT_NEAR(gold_vec[i], meta.sample_vec[i], 1e-12) << i;
        }
    }
}

TEST_F(UniformCtTest, required_parameters_dont_throw)
{
    act::UniformCtData meta;
    amrex::ParmParse pp("Actuator.UniformCtDisk");
    act::utils::ActParser ap("Actuator.UniformCtDisk", "Actuator");
    ASSERT_TRUE(ap.contains("rotor_diameter"));
    ASSERT_NO_THROW(ops::required_parameters(meta, ap));
    EXPECT_DOUBLE_EQ(meta.diameter, 1.0);
    EXPECT_DOUBLE_EQ(meta.epsilon, 1.0);
    EXPECT_DOUBLE_EQ(meta.thrust_coeff[0], 1.0);
    EXPECT_EQ(meta.thrust_coeff.size(), 1);
    EXPECT_EQ(meta.num_force_pts, 3);
}

TEST_F(
    UniformCtTest, single_vel_pos_matches_normal_and_mag_of_specifying_vector)
{
    act::UniformCtData meta;
    amrex::ParmParse pp("Actuator.UniformCtDisk");
    pp.add("yaw", -10.0);
    pp.add("tilt", -10.0);
    pp.add("sample_yaw", 45.0);
    pp.add("sample_tilt", -45.0);
    pp.add("diameters_to_sample", 1.0);
    act::utils::ActParser ap("Actuator.UniformCtDisk", "Actuator");
    ops::required_parameters(meta, ap);
    ops::optional_parameters(meta, ap);

    {
        act::VecList points = {{0, 0, 0}};
        ops::compute_disk_points(
            meta, points, meta.sample_vec, 0, meta.diameters_to_sample);
        amrex::Real mag = vs::mag(points[0]);
        points[0].normalize();
        EXPECT_NEAR(mag, 1.0, 1e-12);
        for (int i = 0; i < 3; i++) {
            EXPECT_NEAR(points[0][i], meta.sample_vec[i], 1e-12);
        }
    }
    {
        act::VecList points = {{0, 0, 0}};
        ops::compute_disk_points(
            meta, points, meta.normal_vec, 0, meta.diameters_to_sample);
        amrex::Real mag = vs::mag(points[0]);
        points[0].normalize();
        EXPECT_NEAR(mag, 1.0, 1e-12);
        for (int i = 0; i < 3; i++) {
            EXPECT_NEAR(points[0][i], meta.normal_vec[i], 1e-12);
        }
    }
}

TEST_F(UniformCtTest, yawed_normal_is_opposite_expected_wind_dir)
{
    // Wind direction
    const vs::Vector north = -vs::Vector::jhat();
    const vs::Vector south = vs::Vector::jhat();
    const vs::Vector east = -vs::Vector::ihat();
    const vs::Vector west = vs::Vector::ihat();

    std::vector<std::pair<vs::Vector, amrex::Real>> couplets;
    couplets.push_back(std::make_pair(north, 0.0));
    couplets.push_back(std::make_pair(south, 180.0));
    couplets.push_back(std::make_pair(west, 270.0));
    couplets.push_back(std::make_pair(east, 90.0));

    amrex::ParmParse pp("Actuator.UniformCtDisk");
    act::utils::ActParser ap("Actuator.UniformCtDisk", "Actuator");
    for (auto&& dir : couplets) {
        act::UniformCtData meta;
        pp.add("yaw", dir.second);
        ops::optional_parameters(meta, ap);
        for (int i = 0; i < AMREX_SPACEDIM; i++) {
            EXPECT_NEAR(dir.first[i], -meta.normal_vec[i], 1e-12)
                << "Failure for yaw: " << dir.second << " index: " << i;
            EXPECT_NEAR(dir.first[i], -meta.sample_vec[i], 1e-12)
                << "Failure for yaw: " << dir.second << " index: " << i;
        }
    }
}

TEST_F(UniformCtTest, sample_yawed_normal_is_opposite_expected_wind_dir)
{
    // Wind direction
    const vs::Vector north = -vs::Vector::jhat();
    const vs::Vector south = vs::Vector::jhat();
    const vs::Vector east = -vs::Vector::ihat();
    const vs::Vector west = vs::Vector::ihat();

    std::vector<std::pair<vs::Vector, amrex::Real>> couplets;
    couplets.push_back(std::make_pair(north, 0.0));
    couplets.push_back(std::make_pair(south, 180.0));
    couplets.push_back(std::make_pair(west, 270.0));
    couplets.push_back(std::make_pair(east, 90.0));

    amrex::ParmParse pp("Actuator.UniformCtDisk");
    act::utils::ActParser ap("Actuator.UniformCtDisk", "Actuator");
    for (auto&& dir : couplets) {
        act::UniformCtData meta;
        pp.add("sample_yaw", dir.second);
        ops::optional_parameters(meta, ap);
        for (int i = 0; i < AMREX_SPACEDIM; i++) {
            EXPECT_NEAR(dir.first[i], -meta.sample_vec[i], 1e-12)
                << "Failure for sample yaw: " << dir.second << " index: " << i;
        }
    }
}
} // namespace amr_wind_tests