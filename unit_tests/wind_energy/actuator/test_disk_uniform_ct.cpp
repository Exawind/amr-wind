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
    pp.add("yaw", 90.0);
    pp.add("tilt", -90.0);
    pp.add("sample_yaw", 45.0);
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

} // namespace amr_wind_tests