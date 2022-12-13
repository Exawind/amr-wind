#ifndef STOKES_WAVES_H_
#define STOKES_WAVES_H_

#include <AMReX_FArrayBox.H>
#include <cmath>

namespace amr_wind::ocean_waves::relaxation_zones {

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void stokes_coefficients(
    int StokesOrder,
    amrex::Real wavenumber,
    amrex::Real waterdepth,
    amrex::Real& c0,
    amrex::Real& a11,
    amrex::Real& a22,
    amrex::Real& b22,
    amrex::Real& c2,
    amrex::Real& d2,
    amrex::Real& e2,
    amrex::Real& a31,
    amrex::Real& a33,
    amrex::Real& b31,
    amrex::Real& a42,
    amrex::Real& a44,
    amrex::Real& b42,
    amrex::Real& b44,
    amrex::Real& c4,
    amrex::Real& d4,
    amrex::Real& e4,
    amrex::Real& a51,
    amrex::Real& a53,
    amrex::Real& a55,
    amrex::Real& b53,
    amrex::Real& b55)
{

    amrex::Real kd = wavenumber * waterdepth;
    if (kd > 50. * M_PI) {
        kd = 50 * M_PI;
    }

    amrex::Real S = 2. * std::exp(2. * kd) / (std::exp(4. * kd) + 1.);
    amrex::Real Sh = std::sinh(kd);
    amrex::Real Th = std::tanh(kd);
    amrex::Real CTh = (1. + std::exp(-2 * kd)) / (1. - std::exp(-2 * kd));

    c0 = std::sqrt(Th);
    a11 = 1. / std::sinh(kd); // Hyperbolic cosecant
    // Second order coefficients
    a22 = 3. * std::pow(S, 2) / (2 * std::pow(1 - S, 2));
    b22 = CTh * (1 + 2. * S) / (2 * (1 - S));
    c2 = std::sqrt(Th) * (2 + 7 * std::pow(S, 2)) / (4 * std::pow(1 - S, 2));
    d2 = -std::sqrt(CTh) / 2.;
    e2 = Th * (2 + 2 * S + 5 * std::pow(S, 2)) / (4 * std::pow(1 - S, 2));
    if (StokesOrder == 2) {
        return;
    }

    // Third order coefficients
    a31 = (-4 - 20 * S + 10 * std::pow(S, 2) - 13 * std::pow(S, 3)) /
          (8 * Sh * std::pow(1 - S, 3));
    a33 = (-2 * std::pow(S, 2) + 11 * std::pow(S, 3)) /
          (8 * Sh * std::pow(1 - S, 3));
    b31 = -3 * (1 + 3 * S + 3 * std::pow(S, 2) + 2 * std::pow(S, 3)) /
          (8 * std::pow(1 - S, 3));
    if (StokesOrder == 3) {
        return;
    }

    // Fourth order coefficients
    a42 = (12 * S - 14 * std::pow(S, 2) - 264 * std::pow(S, 3) -
           45 * std::pow(S, 4) - 13 * std::pow(S, 5)) /
          (24 * std::pow(1 - S, 5));
    a44 = (10 * std::pow(S, 3) - 174 * std::pow(S, 4) + 291 * std::pow(S, 5) +
           278 * std::pow(S, 6)) /
          (48 * (3 + 2 * S) * std::pow(1 - S, 5));
    b42 = CTh *
          (6 - 26 * S - 182 * std::pow(S, 2) - 204 * std::pow(S, 3) -
           25 * std::pow(S, 4) + 26 * std::pow(S, 5)) /
          (6 * (3 + 2 * S) * std::pow(1 - S, 4));
    b44 = CTh *
          (24 + 92 * S + 122 * std::pow(S, 2) + 66 * std::pow(S, 3) +
           67 * std::pow(S, 4) + 34 * std::pow(S, 5)) /
          (24 * (3 + 2 * S) * std::pow(1 - S, 4));
    c4 = std::sqrt(Th) *
         (4 + 32 * S - 116 * std::pow(S, 2) - 400 * std::pow(S, 3) -
          71 * std::pow(S, 4) + 146 * std::pow(S, 5)) /
         (32 * std::pow(1 - S, 5));
    d4 = std::sqrt(CTh) * (2 + 4 * S + std::pow(S, 2) + 2 * std::pow(S, 3)) /
         (8 * std::pow(1 - S, 3));
    e4 = Th *
         (8 + 12 * S - 152 * std::pow(S, 2) - 308 * std::pow(S, 3) -
          42 * std::pow(S, 4) + 77 * std::pow(S, 5)) /
         (32 * std::pow(1 - S, 5));
    if (StokesOrder == 4) {
        return;
    }

    // Fifth order coefficients
    a51 =
        (-1184 + 32 * S + 13232 * std::pow(S, 2) + 21712 * std::pow(S, 3) +
         20940 * std::pow(S, 4) + 12554 * std::pow(S, 5) -
         500 * std::pow(S, 6) - 3341 * std::pow(S, 7) - 670 * std::pow(S, 8)) /
        (64 * Sh * (3 + 2 * S) * (4 + S) * std::pow(1 - S, 6));
    a53 = (4 * S + 105 * pow(S, 2) + 198 * std::pow(S, 3) -
           1376 * std::pow(S, 4) - 1302 * std::pow(S, 5) -
           117 * std::pow(S, 6) + 58 * std::pow(S, 7)) /
          (32 * Sh * (3 + 2 * S) * std::pow(1 - S, 6));
    a55 =
        (-6 * std::pow(S, 3) + 272 * std::pow(S, 4) - 1552 * std::pow(S, 5) +
         852 * std::pow(S, 6) + 2029 * std::pow(S, 7) + 430 * std::pow(S, 8)) /
        (64 * Sh * (3 + 2 * S) * (4 + S) * std::pow(1 - S, 6));
    b53 = 9 *
          (132 + 17 * S - 2216 * std::pow(S, 2) - 5897 * std::pow(S, 3) -
           6292 * std::pow(S, 4) - 2687 * std::pow(S, 5) +
           194 * std::pow(S, 6) + 467 * std::pow(S, 7) + 82 * std::pow(S, 8)) /
          (128 * (3 + 2 * S) * (4 + S) * std::pow(1 - S, 6));
    b55 =
        5 *
        (300 + 1579 * S + 3176 * std::pow(S, 2) + 2949 * std::pow(S, 3) +
         1188 * std::pow(S, 4) + 675 * std::pow(S, 5) + 1326 * std::pow(S, 6) +
         827 * std::pow(S, 7) + 130 * std::pow(S, 8)) /
        (384 * (3 + 2 * S) * (4 + S) * std::pow(1 - S, 6));
    if (StokesOrder == 5) {
        return;
    }

    if (StokesOrder > 5 || StokesOrder < 2) {
        amrex::Abort(
            "invalid stokes order speficied. It should be between 2,3,4 or 5 ");
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real my_sinh_sin(
    int i,
    int j,
    amrex::Real phase,
    amrex::Real a11,
    amrex::Real a22,
    amrex::Real a31,
    amrex::Real a33,
    amrex::Real a42,
    amrex::Real a44,
    amrex::Real a51,
    amrex::Real a53,
    amrex::Real a55,
    amrex::Real eps,
    amrex::Real k,
    amrex::Real waterdepth)
{
    amrex::Real D = 0.0;
    if (i == 1 && j == 1) {
        D = a11;
    }
    if (i == 2 && j == 2) {
        D = a22;
    }
    if (i == 3 && j == 1) {
        D = a31;
    }
    if (i == 3 && j == 3) {
        D = a33;
    }
    if (i == 4 && j == 2) {
        D = a42;
    }
    if (i == 4 && j == 4) {
        D = a44;
    }
    if (i == 5 && j == 1) {
        D = a51;
    }
    if (i == 5 && j == 3) {
        D = a53;
    }
    if (i == 5 && j == 5) {
        D = a55;
    }

    return std::pow(eps, i) * D * j * k * std::sinh(j * k * waterdepth) *
           std::sin(j * phase);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real my_cosh_cos(
    int i,
    int j,
    amrex::Real phase,
    amrex::Real a11,
    amrex::Real a22,
    amrex::Real a31,
    amrex::Real a33,
    amrex::Real a42,
    amrex::Real a44,
    amrex::Real a51,
    amrex::Real a53,
    amrex::Real a55,
    amrex::Real eps,
    amrex::Real k,
    amrex::Real waterdepth)
{
    amrex::Real D = 0.0;
    if (i == 1 && j == 1) {
        D = a11;
    }
    if (i == 2 && j == 2) {
        D = a22;
    }
    if (i == 3 && j == 1) {
        D = a31;
    }
    if (i == 3 && j == 3) {
        D = a33;
    }
    if (i == 4 && j == 2) {
        D = a42;
    }
    if (i == 4 && j == 4) {
        D = a44;
    }
    if (i == 5 && j == 1) {
        D = a51;
    }
    if (i == 5 && j == 3) {
        D = a53;
    }
    if (i == 5 && j == 5) {
        D = a55;
    }

    return std::pow(eps, i) * D * j * k * std::cosh(j * k * waterdepth) *
           std::cos(j * phase);
}

// Based on Fenton 1985
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void stokes_waves(
    int StokesOrder,
    amrex::Real wavelength,
    amrex::Real waterdepth,
    amrex::Real waveheight,
    amrex::Real x,
    amrex::Real time,
    amrex::Real& eta,
    amrex::Real& u_w,
    amrex::Real& v_w,
    amrex::Real& w_w)
{
    const amrex::Real k = 2.0 * M_PI / wavelength;

    // some parameters
    amrex::Real c0{0.0};
    amrex::Real a11{0.0}, a22{0.0}, b22{0.0}, c2{0.0}, d2{0.0}, e2{0.0};
    amrex::Real a31{0.0}, a33{0.0}, b31{0.0}, a42{0.0}, a44{0.0};
    amrex::Real b42{0.0}, b44{0.0}, c4{0.0}, d4{0.0}, e4{0.0};
    amrex::Real a51{0.0}, a53{0.0}, a55{0.0}, b53{0.0}, b55{0.0};

    stokes_coefficients(
        StokesOrder, k, waterdepth, c0, a11, a22, b22, c2, d2, e2, a31, a33,
        b31, a42, a44, b42, b44, c4, d4, e4, a51, b53, a55, b53, b55);

    const amrex::Real eps = k * waveheight / 2.; // Steepness (ka)
    const amrex::Real c = (c0 + std::pow(eps, 2) * c2 + std::pow(eps, 4) * c4) *
                          std::sqrt(9.81 / k);
    // const amrex::Real Q =
    //    c * waterdepth * std::sqrt(std::pow(k, 3) / 9.81) +
    //    d2 * std::pow(eps, 2) +
    //    d4 * std::pow(eps, 4) * std::sqrt(9.81 / std::pow(k, 3));

    const amrex::Real omega = c * k;
    const amrex::Real phase = k * x - omega * time;

    eta =
        (eps * std::cos(phase)                           // first order term
         + std::pow(eps, 2) * b22 * std::cos(2. * phase) // second order term
         + std::pow(eps, 3) * b31 * (std::cos(phase) - std::cos(3. * phase)) +
         std::pow(eps, 4) * b42 *
             (std::cos(2. * phase) + b44 * std::cos(4. * phase)) +
         std::pow(eps, 5) *
             (-(b53 + b55) * std::cos(phase) + b53 * std::cos(3 * phase) +
              b55 * std::cos(5 * phase)));

    // Compute cosines
    const amrex::Real cc11 = my_cosh_cos(
        1, 1, phase, a11, a22, a31, a33, a42, a44, a51, a53, a55, eps, k,
        waterdepth);
    const amrex::Real cc22 = my_cosh_cos(
        2, 2, phase, a11, a22, a31, a33, a42, a44, a51, a53, a55, eps, k,
        waterdepth);
    const amrex::Real cc31 = my_cosh_cos(
        3, 1, phase, a11, a22, a31, a33, a42, a44, a51, a53, a55, eps, k,
        waterdepth);
    const amrex::Real cc33 = my_cosh_cos(
        3, 3, phase, a11, a22, a31, a33, a42, a44, a51, a53, a55, eps, k,
        waterdepth);
    const amrex::Real cc42 = my_cosh_cos(
        4, 2, phase, a11, a22, a31, a33, a42, a44, a51, a53, a55, eps, k,
        waterdepth);
    const amrex::Real cc44 = my_cosh_cos(
        4, 4, phase, a11, a22, a31, a33, a42, a44, a51, a53, a55, eps, k,
        waterdepth);
    const amrex::Real cc51 = my_cosh_cos(
        5, 1, phase, a11, a22, a31, a33, a42, a44, a51, a53, a55, eps, k,
        waterdepth);
    const amrex::Real cc53 = my_cosh_cos(
        5, 3, phase, a11, a22, a31, a33, a42, a44, a51, a53, a55, eps, k,
        waterdepth);
    const amrex::Real cc55 = my_cosh_cos(
        5, 5, phase, a11, a22, a31, a33, a42, a44, a51, a53, a55, eps, k,
        waterdepth);

    // Compute sines
    const amrex::Real ss11 = my_sinh_sin(
        1, 1, phase, a11, a22, a31, a33, a42, a44, a51, a53, a55, eps, k,
        waterdepth);
    const amrex::Real ss22 = my_sinh_sin(
        2, 2, phase, a11, a22, a31, a33, a42, a44, a51, a53, a55, eps, k,
        waterdepth);
    const amrex::Real ss31 = my_sinh_sin(
        3, 1, phase, a11, a22, a31, a33, a42, a44, a51, a53, a55, eps, k,
        waterdepth);
    const amrex::Real ss33 = my_sinh_sin(
        3, 3, phase, a11, a22, a31, a33, a42, a44, a51, a53, a55, eps, k,
        waterdepth);
    const amrex::Real ss42 = my_sinh_sin(
        4, 2, phase, a11, a22, a31, a33, a42, a44, a51, a53, a55, eps, k,
        waterdepth);
    const amrex::Real ss44 = my_sinh_sin(
        4, 4, phase, a11, a22, a31, a33, a42, a44, a51, a53, a55, eps, k,
        waterdepth);
    const amrex::Real ss51 = my_sinh_sin(
        5, 1, phase, a11, a22, a31, a33, a42, a44, a51, a53, a55, eps, k,
        waterdepth);
    const amrex::Real ss53 = my_sinh_sin(
        5, 3, phase, a11, a22, a31, a33, a42, a44, a51, a53, a55, eps, k,
        waterdepth);
    const amrex::Real ss55 = my_sinh_sin(
        5, 5, phase, a11, a22, a31, a33, a42, a44, a51, a53, a55, eps, k,
        waterdepth);

    u_w = c0 * std::sqrt(9.81 / std::pow(k, 3)) *
          (cc11 + cc22 + cc31 + cc33 + cc42 + cc44 + cc51 + cc53 + cc55);
    v_w = 0.0;
    w_w = c0 * std::sqrt(9.81 / std::pow(k, 3)) *
          (ss11 + ss22 + ss31 + ss33 + ss42 + ss44 + ss51 + ss53 + ss55);
}

} // namespace amr_wind::ocean_waves::relaxation_zones
#endif
