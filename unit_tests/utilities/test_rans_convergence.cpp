#include <cmath>
#include <limits>
#include <vector>

#include "ks_test_utils/AmrexTest.H"
#include "src/utilities/output_quantities/RANSConvergence.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace kynema_sgf_tests {

namespace {

namespace rc = kynema_sgf::rans_convergence;

//! Build an exponentially decaying envelope history
std::vector<amrex::Real> decaying_series(
    const std::vector<amrex::Real>& times,
    const amrex::Real amplitude,
    const amrex::Real rate)
{
    std::vector<amrex::Real> vals;
    vals.reserve(times.size());
    for (const auto t : times) {
        vals.push_back(amplitude * std::exp(-rate * t));
    }
    return vals;
}

std::vector<amrex::Real>
uniform_times(const int n, const amrex::Real t0, const amrex::Real dt)
{
    std::vector<amrex::Real> times;
    times.reserve(n);
    for (int i = 0; i < n; ++i) {
        times.push_back(t0 + (i * dt));
    }
    return times;
}

//! Absolute tolerance for comparing a value against an exact expectation.
//!
//! Scaled to the working precision so that these tests mean the same thing in
//! single and double precision. The fit takes logarithms and an exponential,
//! so a few thousand epsilon of relative error is expected rather than a few.
amrex::Real close(const amrex::Real expected)
{
    constexpr amrex::Real factor = 1.0e4_rt;
    return factor * std::numeric_limits<amrex::Real>::epsilon() *
           std::max(std::abs(expected), 1.0_rt);
}

} // namespace

TEST(RANSConvergence, effective_tolerance_takes_the_larger_term)
{
    // Relative term dominates for a large mean
    EXPECT_NEAR(
        rc::RANSConvergence::effective_tolerance(10.0_rt, 0.01_rt, 0.01_rt),
        0.1_rt, close(0.1_rt));
    // Absolute floor takes over as the mean approaches zero
    EXPECT_NEAR(
        rc::RANSConvergence::effective_tolerance(0.0_rt, 0.01_rt, 0.01_rt),
        0.01_rt, close(0.01_rt));
    // A negative mean is treated by magnitude
    EXPECT_NEAR(
        rc::RANSConvergence::effective_tolerance(-10.0_rt, 0.01_rt, 0.01_rt),
        0.1_rt, close(0.1_rt));
}

TEST(RANSConvergence, envelope_fit_recovers_a_known_decay)
{
    const amrex::Real amplitude = 5.0_rt;
    const amrex::Real rate = 1.0e-3_rt;
    // Kept short enough that the series is still above the threshold at the
    // last sample, which is the only regime in which extrapolating forward
    // means anything
    const auto times = uniform_times(40, 0.0_rt, 10.0_rt);
    const auto vals = decaying_series(times, amplitude, rate);

    const auto fit =
        rc::RANSConvergence::fit_envelope_decay(times, vals, 1.0_rt, 5);

    ASSERT_TRUE(fit.valid);
    EXPECT_NEAR(fit.rate, rate, close(rate));
    EXPECT_NEAR(fit.amplitude, amplitude, close(amplitude));
    EXPECT_NEAR(fit.rsq, 1.0_rt, close(1.0_rt));

    // A exp(-rate t) = 1 gives t = ln(A)/rate, measured from the last sample
    const amrex::Real expected = (std::log(amplitude) / rate) - times.back();
    EXPECT_NEAR(fit.time_to_threshold, expected, close(expected));
}

TEST(RANSConvergence, envelope_fit_rejects_a_growing_envelope)
{
    const auto times = uniform_times(20, 0.0_rt, 100.0_rt);
    // Negative rate means the spread is growing, so it never reaches the
    // threshold from above
    const auto vals = decaying_series(times, 0.5_rt, -1.0e-3_rt);

    const auto fit =
        rc::RANSConvergence::fit_envelope_decay(times, vals, 1.0_rt, 5);

    EXPECT_FALSE(fit.valid);
    EXPECT_LT(fit.rate, 0.0_rt);
}

TEST(RANSConvergence, envelope_fit_rejects_too_few_samples)
{
    const auto times = uniform_times(4, 0.0_rt, 100.0_rt);
    const auto vals = decaying_series(times, 5.0_rt, 1.0e-3_rt);

    const auto fit =
        rc::RANSConvergence::fit_envelope_decay(times, vals, 1.0_rt, 5);

    EXPECT_FALSE(fit.valid);
}

TEST(RANSConvergence, envelope_fit_rejects_non_positive_spreads)
{
    const auto times = uniform_times(10, 0.0_rt, 100.0_rt);
    auto vals = decaying_series(times, 5.0_rt, 1.0e-3_rt);
    // An identical pair of samples collapses the spread to exactly zero,
    // which has no logarithm
    vals[4] = 0.0_rt;

    const auto fit =
        rc::RANSConvergence::fit_envelope_decay(times, vals, 1.0_rt, 5);

    EXPECT_FALSE(fit.valid);
}

TEST(RANSConvergence, envelope_fit_rejects_a_threshold_already_met)
{
    const auto times = uniform_times(20, 0.0_rt, 100.0_rt);
    // The series decays from 0.5, so it is already below a threshold of one
    const auto vals = decaying_series(times, 0.5_rt, 1.0e-3_rt);

    const auto fit =
        rc::RANSConvergence::fit_envelope_decay(times, vals, 1.0_rt, 5);

    EXPECT_FALSE(fit.valid);
    EXPECT_LT(fit.time_to_threshold, 0.0_rt);
}

TEST(RANSConvergence, envelope_fit_reports_a_poor_fit_through_rsq)
{
    // A spread that has stalled on a noise floor is not an exponential decay.
    // The fit still returns, but with an r-squared low enough that the caller
    // will discard the estimate rather than quote a confident wrong number.
    const auto times = uniform_times(20, 0.0_rt, 100.0_rt);
    std::vector<amrex::Real> vals(times.size(), 2.0_rt);
    for (size_t i = 0; i < vals.size(); ++i) {
        vals[i] += ((i % 2 == 0) ? 0.1_rt : -0.1_rt);
    }

    const auto fit =
        rc::RANSConvergence::fit_envelope_decay(times, vals, 1.0_rt, 5);

    EXPECT_LT(fit.rsq, 0.5_rt);
}

TEST(RANSConvergence, envelope_fit_tolerates_noise_on_a_real_decay)
{
    // The estimate does not need to be precise, only usable. A ten percent
    // multiplicative wobble should still recover the decay rate closely
    // enough to be worth printing.
    const amrex::Real amplitude = 8.0_rt;
    const amrex::Real rate = 5.0e-4_rt;
    const auto times = uniform_times(60, 0.0_rt, 20.0_rt);
    auto vals = decaying_series(times, amplitude, rate);
    for (size_t i = 0; i < vals.size(); ++i) {
        vals[i] *= 1.0_rt + ((i % 3 == 0) ? 0.1_rt : -0.05_rt);
    }

    const auto fit =
        rc::RANSConvergence::fit_envelope_decay(times, vals, 1.0_rt, 5);

    ASSERT_TRUE(fit.valid);
    EXPECT_NEAR(fit.rate, rate, 0.1_rt * rate);
    // Noise of this size costs some of the fit quality, but the result stays
    // well clear of the default eta_min_rsq gate of 0.5, so the estimate
    // would be reported rather than discarded
    EXPECT_GT(fit.rsq, 0.8_rt);
}

} // namespace kynema_sgf_tests
