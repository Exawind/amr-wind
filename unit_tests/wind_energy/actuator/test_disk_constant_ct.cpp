#include "aw_test_utils/AmrexTest.H"
#include "test_act_utils.H"

#include "amr-wind/wind_energy/actuator/disk/constant_ct_ops.H"

namespace amr_wind_tests {

namespace {
class ConstantCtTest : public AmrexTest
{
    void SetUp() override
    {
        AmrexTest::SetUp();
        amrex::ParmParse pp("ConstantCtTest");
    }
};
} // namespace

TEST_F(ConstantCtTest, dummy_test) { ASSERT_TRUE(false); }

} // namespace amr_wind_tests