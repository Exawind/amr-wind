#include "aw_test_utils/AmrexTest.H"

#include "amr-wind/utilities/linear_interpolation.H"
#include "amr-wind/utilities/trig_ops.H"
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

TEST(LinearInterpolation, nearest_search)
{
    namespace interp = amr_wind::interp;
    std::vector<amrex::Real> xvec(10);
    std::iota(xvec.begin(), xvec.end(), 0.0);

    const auto* start = xvec.data();
    const auto* end = xvec.data() + xvec.size();
    {
        const auto idx = interp::nearest_search(start, end, 5.0);
        EXPECT_EQ(idx.idx, 5);
        EXPECT_EQ(idx.lim, interp::Limits::VALID);
    }
    {
        const auto idx = interp::nearest_search(start, end, 1.1);
        EXPECT_EQ(idx.idx, 1);
        EXPECT_EQ(idx.lim, interp::Limits::VALID);
    }
    {
        const auto idx = interp::nearest_search(start, end, 1.6);
        EXPECT_EQ(idx.idx, 2);
        EXPECT_EQ(idx.lim, interp::Limits::VALID);
    }
    {
        const auto idx = interp::nearest_search(start, end, 3.5);
        EXPECT_EQ(idx.idx, 3);
        EXPECT_EQ(idx.lim, interp::Limits::VALID);
    }
    {
        const auto idx = interp::nearest_search(start, end, -1.0);
        EXPECT_EQ(idx.idx, 0);
        EXPECT_EQ(idx.lim, interp::Limits::LOWLIM);
    }
    {
        const auto idx = interp::nearest_search(start, end, 9.1);
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

TEST(LinearInterpolation, lin_interp_single_multicomponent)
{
    namespace interp = amr_wind::interp;

    const int ncomp = 3;
    std::vector<amrex::Real> xvec(10), yvec(ncomp * 10);
    std::iota(xvec.begin(), xvec.end(), 0.0);
    const amrex::Vector<amrex::Real> mult_facs = {
        2.0 + 10.0 * amrex::Random(), 2.0 + 10.0 * amrex::Random(),
        2.0 + 10.0 * amrex::Random()};
    for (int i = 0; i < static_cast<int>(xvec.size()); i++) {
        for (int n = 0; n < ncomp; n++) {
            yvec[ncomp * i + n] = mult_facs[n] * xvec[i];
        }
    }

    std::vector<amrex::Real> xtest{2.5, 4.5, 6.3, 8.8};
    for (const auto& x : xtest) {
        for (int n = 0; n < ncomp; n++) {
            const auto y = interp::linear(xvec, yvec, x, 3, n);
            EXPECT_NEAR(y, mult_facs[n] * x, 1.0e-12);
        }
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

TEST(LinearInterpolation, lin_interp_angle)
{
    namespace interp = amr_wind::interp;

    const int vecsize = 8;
    std::vector<amrex::Real> xvec(vecsize);
    std::iota(xvec.begin(), xvec.end(), 0.0);
    // Set up with angles that are awkward to interpolate between
    std::vector<amrex::Real> yvec_deg{0.0,   30.0,   330.0, 300.0,
                                      -15.0, -180.0, 120.0, -150.0};
    // Create duplicate list in radians
    std::vector<amrex::Real> yvec_rad(vecsize);
    for (size_t i = 0; i < yvec_rad.size(); ++i) {
        yvec_rad[i] = amr_wind::utils::radians(yvec_deg[i]);
    }

    // Set up output vectors, interp locations, and golds
    std::vector<amrex::Real> xtest(vecsize), ytest_deg(vecsize),
        ytest_rad(vecsize);
    std::iota(xtest.begin(), xtest.end(), 0.5);
    xtest[7] = 4.0;
    std::vector<amrex::Real> ygold_deg{15.,   0.,   315., 322.5,
                                       -97.5, -210, 165., 360. - 15.};

    // Upper bound is 360 because y is in degrees
    interp::linear_angle(xvec, yvec_deg, xtest, ytest_deg, 360.);
    // Upper bound is 2pi because y is in radians
    interp::linear_angle(xvec, yvec_rad, xtest, ytest_rad, 2. * M_PI);
    for (size_t i = 0; i < xtest.size(); ++i) {
        EXPECT_NEAR(ytest_deg[i], ygold_deg[i], 1.0e-12);
        EXPECT_NEAR(
            ytest_rad[i], amr_wind::utils::radians(ygold_deg[i]), 1.0e-12);
    }
}

} // namespace amr_wind_tests
