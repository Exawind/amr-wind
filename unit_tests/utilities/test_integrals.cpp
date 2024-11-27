#include "aw_test_utils/AmrexTest.H"
#include "amr-wind/utilities/integrals.H"
#include "amr-wind/utilities/constants.H"

namespace amr_wind_tests {

TEST(Integrals, trapezoid_integration)
{
    const amrex::Real xa = 0.0;
    const amrex::Real xb = 1.2;
    const int n = 10000;

    const amrex::Real integ = amr_wind::utils::trapz(
        xa, xb, n,
        [] AMREX_GPU_DEVICE(const amrex::Real x) noexcept { return x * x; });
    EXPECT_NEAR(integ, 0.576, amr_wind::constants::LOOSE_TOL);
}

} // namespace amr_wind_tests
