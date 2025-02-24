#include <AMReX_Gpu.H>
#include "aw_test_utils/AmrexTest.H"
#include "amr-wind/eos_models/EOSModel.H"

namespace amr_wind_tests {

void test_eos_impl()
{
    const auto eos = amr_wind::eos::GammaLaw(1.01325e5);
    amrex::Gpu::DeviceScalar<amrex::Real> val(0.0);
    auto* d_val = val.dataPtr();

    amrex::ParallelFor(1, [=] AMREX_GPU_DEVICE(int /*unused*/) {
        d_val[0] = eos.p_rth(1.225, 300.0, 0.5);
    });
    EXPECT_NEAR(
        val.dataValue(), 244859.65251771925, amr_wind::constants::LOOSE_TOL);

    amrex::ParallelFor(1, [=] AMREX_GPU_DEVICE(int /*unused*/) {
        d_val[0] = eos.dp_constanttheta(1.225, 300.0, 0.5);
    });
    EXPECT_NEAR(
        val.dataValue(), 279839.60287739342, amr_wind::constants::LOOSE_TOL);
}

TEST(EOSModelTest, test_eos) { test_eos_impl(); }
} // namespace amr_wind_tests
