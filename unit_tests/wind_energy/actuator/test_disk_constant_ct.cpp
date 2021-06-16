#include "aw_test_utils/AmrexTest.H"
#include "test_act_utils.H"

#include "amr-wind/wind_energy/actuator/disk/constant_ct_ops.H"

namespace amr_wind_tests {

namespace {
class ConstantCtTest : public AmrexTest
{
protected:
    void SetUp() override
    {
        AmrexTest::SetUp();
        amrex::ParmParse pp("Actuator");
    }
};
} // namespace

namespace act = amr_wind::actuator;
namespace ops = act::ops;
namespace vs = amr_wind::vs;

TEST_F(ConstantCtTest, compute_vecs_from_yaw)
{
    act::ConstantCtData meta;
    amrex::ParmParse pp("Actuator.ConstantCtDisk");
    pp.add("yaw", 90.0);
    pp.add("sample_yaw", 45.0);
    act::utils::ActParser ap("Actuator.ConstantCtDisk", "Actuator");
    ASSERT_TRUE(ap.contains("yaw"));
    ops::optional_parameters(meta, ap);
    {
        const vs::Vector gold_nom = vs::Vector::jhat();
        for (int i = 0; i < 3; i++) {
            EXPECT_NEAR(gold_nom[i], meta.normal_vec[i], 1e-12);
        }
    }
    {
        vs::Vector gold_nom = vs::Vector::jhat() + vs::Vector::ihat();
        gold_nom.normalize();

        for (int i = 0; i < 3; i++) {
            EXPECT_NEAR(gold_nom[i], meta.sample_vec[i], 1e-12);
        }
    }
}

} // namespace amr_wind_tests