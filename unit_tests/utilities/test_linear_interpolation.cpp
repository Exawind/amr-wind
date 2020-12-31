#include "aw_test_utils/AmrexTest.H"

#include "amr-wind/utilities/linear_interpolation.H"
#include "AMReX_Random.H"

#include <vector>
#include <numeric>

namespace amr_wind_tests {

TEST(LinearInterpolation, check_bounds)
{
    namespace interp = amr_wind::interp;
    std::vector<amrex::Real> xvec(10);
    std::iota(xvec.begin(), xvec.end(), 0.0);

    const auto* start = xvec.data();
    const auto* end = xvec.data() + xvec.size();
    {
        const auto idx = interp::check_bounds(start, end, -1.0);
        EXPECT_EQ(idx.idx, 0);
        EXPECT_EQ(idx.lim, interp::Limits::LOWLIM);
    }
    {
        const auto idx = interp::check_bounds(start, end, 9.1);
        EXPECT_EQ(idx.idx, 9);
        EXPECT_EQ(idx.lim, interp::Limits::UPLIM);
    }
    {
        const amrex::Real xinp = 9.0 * amrex::Random();
        const auto idx = interp::check_bounds(start, end, xinp);
        EXPECT_EQ(idx.idx, 0);
        EXPECT_EQ(idx.lim, interp::Limits::VALID);
    }
}

TEST(LinearInterpolation, bisection_search)
{
    namespace interp = amr_wind::interp;
    std::vector<amrex::Real> xvec(10);
    std::iota(xvec.begin(), xvec.end(), 0.0);

    const auto* start = xvec.data();
    const auto* end = xvec.data() + xvec.size();
    {
        const auto idx = interp::bisection_search(start, end, 5.0);
        EXPECT_EQ(idx.idx, 4);
        EXPECT_EQ(idx.lim, interp::Limits::VALID);
    }
    {
        const auto idx = interp::bisection_search(start, end, -1.0);
        EXPECT_EQ(idx.idx, 0);
        EXPECT_EQ(idx.lim, interp::Limits::LOWLIM);
    }
    {
        const auto idx = interp::bisection_search(start, end, 9.1);
        EXPECT_EQ(idx.idx, 9);
        EXPECT_EQ(idx.lim, interp::Limits::UPLIM);
    }
}

TEST(LinearInterpolation, find_index)
{
    namespace interp = amr_wind::interp;
    std::vector<amrex::Real> xvec(10);
    std::iota(xvec.begin(), xvec.end(), 0.0);

    const auto* start = xvec.data();
    const auto* end = xvec.data() + xvec.size();
    {
        const auto idx = interp::find_index(start, end, 5.0);
        EXPECT_EQ(idx.idx, 4);
        EXPECT_EQ(idx.lim, interp::Limits::VALID);
    }
    {
        const auto idx = interp::find_index(start, end, 5.0, 3);
        EXPECT_EQ(idx.idx, 4);
        EXPECT_EQ(idx.lim, interp::Limits::VALID);
    }
    {
        const auto idx = interp::find_index(start, end, -1.0);
        EXPECT_EQ(idx.idx, 0);
        EXPECT_EQ(idx.lim, interp::Limits::LOWLIM);
    }
    {
        const auto idx = interp::find_index(start, end, 9.1);
        EXPECT_EQ(idx.idx, 9);
        EXPECT_EQ(idx.lim, interp::Limits::UPLIM);
    }
}

TEST(LinearInterpolation, lin_interp_single)
{
    namespace interp = amr_wind::interp;

    const amrex::Real mult_fac = 2.0 + 10.0 * amrex::Random();
    std::vector<amrex::Real> xvec(10), yvec(10);
    std::iota(xvec.begin(), xvec.end(), 0.0);
    std::transform(
        xvec.begin(), xvec.end(), yvec.begin(),
        [mult_fac](const amrex::Real x) { return mult_fac * x; });

    std::vector<amrex::Real> xtest{2.5, 4.5, 6.3, 8.8};
    for (const auto& x : xtest) {
        const auto y = interp::linear(xvec, yvec, x);
        EXPECT_NEAR(y, mult_fac * x, 1.0e-12);
    }
}

TEST(LinearInterpolation, lin_interp)
{
    namespace interp = amr_wind::interp;

    const amrex::Real mult_fac = 2.0 + 10.0 * amrex::Random();
    std::vector<amrex::Real> xvec(10), yvec(10);
    std::iota(xvec.begin(), xvec.end(), 0.0);
    std::transform(
        xvec.begin(), xvec.end(), yvec.begin(),
        [mult_fac](const amrex::Real x) { return mult_fac * x; });

    std::vector<amrex::Real> xtest{2.5, 4.5, 6.3, 8.8};
    std::vector<amrex::Real> ytest(xtest.size());

    interp::linear(xvec, yvec, xtest, ytest);
    for (size_t i = 0; i < xtest.size(); ++i) {
        EXPECT_NEAR(ytest[i], mult_fac * xtest[i], 1.0e-12);
    }
}

TEST(LinearInterpolation, lin_monotonic)
{
    namespace interp = amr_wind::interp;

    const amrex::Real mult_fac = 2.0 + 10.0 * amrex::Random();
    std::vector<amrex::Real> xvec(10), yvec(10);
    std::iota(xvec.begin(), xvec.end(), 0.0);
    std::transform(
        xvec.begin(), xvec.end(), yvec.begin(),
        [mult_fac](const amrex::Real x) { return mult_fac * x; });

    std::vector<amrex::Real> xtest{2.5, 4.5, 6.3, 8.8};
    std::vector<amrex::Real> ytest(xtest.size());

    interp::linear_monotonic(xvec, yvec, xtest, ytest);
    for (size_t i = 0; i < xtest.size(); ++i) {
        EXPECT_NEAR(ytest[i], mult_fac * xtest[i], 1.0e-12);
    }
}

} // namespace amr_wind_tests
