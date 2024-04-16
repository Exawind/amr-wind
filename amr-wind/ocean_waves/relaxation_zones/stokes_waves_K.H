#ifndef STOKES_WAVES_H_
#define STOKES_WAVES_H_

#include <AMReX_FArrayBox.H>
#include <cmath>

// Stokes waves theory adapted from
// Fenton, J., Fifth Order Stokes Theory for Steady Waves
// Journal of Waterway, Port, Coastal and Ocean Engineering, 1985, 111, 216-234

namespace amr_wind::ocean_waves::relaxation_zones {

// Compute wavelength as a function of wave period, water depth, and g
AMREX_FORCE_INLINE amrex::Real stokes_wave_length(
    const amrex::Real T,
    const amrex::Real d,
    const amrex::Real H,
    const int order,
    const amrex::Real g,
    const amrex::Real tol,
    const int iter_max)
{
    // Calculate constants and derivatives that do not change with iteration
    const amrex::Real omega = 2.0 * M_PI / T;
    const amrex::Real depsdk = H / 2.0;

    // Begin Newton-Raphson iterations
    int iter = 0;
    // First guess is first-order
    amrex::Real k = omega / std::sqrt(g);
    // Cannot skip loop
    amrex::Real f = tol + 1.0;
    while (std::abs(f) > tol && iter < iter_max) {
        // Calculate current constants
        const amrex::Real S = 1.0 / std::cosh(2.0 * k * d);
        const amrex::Real C = 1.0 - S;
        const amrex::Real eps = k * H / 2.0;
        const amrex::Real C0 = std::sqrt(std::tanh(k * d));
        const amrex::Real C2 = C0 * (2.0 + 7.0 * S * S) / (4.0 * C * C);
        const amrex::Real numterm_C4 =
            (4.0 + 32.0 * S - 116.0 * std::pow(S, 2) - 400.0 * std::pow(S, 3) -
             71.0 * std::pow(S, 4) + 146.0 * std::pow(S, 5));
        const amrex::Real C4 = C0 * numterm_C4 / (32.0 * std::pow(C, 5));
        // Calculate pure derivates
        const amrex::Real dSdk = -2.0 * d * std::sinh(2.0 * k * d) /
                                 std::pow(std::cosh(2.0 * k * d), 2);
        const amrex::Real dCdk = -dSdk;
        const amrex::Real dC0dk =
            d / (2.0 * C0 * std::pow(std::cosh(k * d), 2));
        // Calculate derivatives with products
        const amrex::Real dC2dk =
            (4 * std::pow(C, 2) *
                 (dC0dk * (2.0 + 7.0 * std::pow(S, 2)) + C0 * 14.0 * S * dSdk) -
             C0 * (2.0 + 7.0 * std::pow(S, 2)) * 8.0 * C * dCdk) /
            (16.0 * std::pow(C, 4));
        const amrex::Real dC4dk =
            (32.0 * std::pow(C, 5) *
                 (dC0dk * numterm_C4 +
                  C0 * (32.0 * dSdk - 232.0 * S * dSdk -
                        1200 * std::pow(S, 2) * dSdk -
                        284.0 * std::pow(S, 3) * dSdk + 730 * std::pow(S, 4))) -
             C0 * numterm_C4 * 160.0 * std::pow(C, 4) * dCdk) /
            (1024.0 * std::pow(C, 10));

        // Calculate derivative for loop convergence
        amrex::Real dfdk = g * std::pow(C0, 2) + g * k * 2.0 * C0 * dC0dk;
        // Add additional terms depending on order
        if (order >= 2) {
            dfdk +=
                g * (2.0 * C0 * std::pow(eps, 2) * C2 +
                     std::pow(eps, 4) * std::pow(C2, 2)) +
                g * k *
                    (2.0 * dC0dk * std::pow(eps, 2) * C2 +
                     2.0 * C0 *
                         (2.0 * eps * depsdk * C2 + std::pow(eps, 2) * dC2dk) +
                     4.0 * std::pow(eps, 3) * std::pow(C2, 2) * depsdk +
                     std::pow(eps, 4) * 2.0 * C2 * dC2dk);
        }
        if (order >= 4) {
            dfdk += g * (2.0 * std::pow(eps, 4) * C0 * C4 +
                         2.0 * std::pow(eps, 6) * C2 * C4 +
                         std::pow(eps, 8) * std::pow(C4, 2)) +
                    g * k *
                        (8.0 * std::pow(eps, 3) * depsdk * C0 * C4 +
                         2.0 * std::pow(eps, 4) * (C0 * dC4dk + C4 * dC0dk) +
                         12 * std::pow(eps, 5) * depsdk * C2 * C4 +
                         2.0 * std::pow(eps, 6) * (C2 * dC4dk + C4 * dC2dk) +
                         8.0 * std::pow(eps, 7) * depsdk * std::pow(C4, 2) +
                         std::pow(eps, 8) * 2.0 * C4 * dC4dk);
        }
        k = k - f / dfdk;

        iter += 1;
        f = g * k * std::pow(C0, 2);
        // Add additional terms depending on order
        if (order >= 2) {
            f += g * k *
                 (2.0 * C0 * std::pow(eps, 2) * C2 +
                  std::pow(eps, 4) * std::pow(C2, 2));
        }
        if (order >= 4) {
            f += g * k *
                 (2.0 * std::pow(eps, 4) * C0 * C4 +
                  2.0 * std::pow(eps, 6) * C2 * C4 +
                  std::pow(eps, 8) * std::pow(C4, 2));
        }
        // Subtract omega^2 to measure convergence
        f -= omega * omega;
    }

    if (k < tol) {
        // Return negative wavelength if faulty result
        return -1;
    }
    if (std::isnan(k)) {
        return -2;
    }
    if (iter == iter_max) {
        return -3;
    }
    // Return wavelength calculated from wavenumber
    return 2.0 * M_PI / k;
}

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
            "invalid stokes order specified. It should be between 2,3,4 or 5 ");
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real stokes_sinh_sin(
    int m,
    int n,
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
    amrex::Real wavenumber,
    amrex::Real waterdepth,
    amrex::Real zsl,
    amrex::Real z)
{
    amrex::Real D = 0.0;
    if (m == 1 && n == 1) {
        D = a11;
    }
    if (m == 2 && n == 2) {
        D = a22;
    }
    if (m == 3 && n == 1) {
        D = a31;
    }
    if (m == 3 && n == 3) {
        D = a33;
    }
    if (m == 4 && n == 2) {
        D = a42;
    }
    if (m == 4 && n == 4) {
        D = a44;
    }
    if (m == 5 && n == 1) {
        D = a51;
    }
    if (m == 5 && n == 3) {
        D = a53;
    }
    if (m == 5 && n == 5) {
        D = a55;
    }

    return std::pow(eps, m) * D * n * wavenumber *
           std::sinh(n * wavenumber * (waterdepth + (z - zsl))) *
           std::sin(n * phase);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real stokes_cosh_cos(
    int m,
    int n,
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
    amrex::Real wavenumber,
    amrex::Real waterdepth,
    amrex::Real zsl,
    amrex::Real z)
{
    amrex::Real D = 0.0;
    if (m == 1 && n == 1) {
        D = a11;
    }
    if (m == 2 && n == 2) {
        D = a22;
    }
    if (m == 3 && n == 1) {
        D = a31;
    }
    if (m == 3 && n == 3) {
        D = a33;
    }
    if (m == 4 && n == 2) {
        D = a42;
    }
    if (m == 4 && n == 4) {
        D = a44;
    }
    if (m == 5 && n == 1) {
        D = a51;
    }
    if (m == 5 && n == 3) {
        D = a53;
    }
    if (m == 5 && n == 5) {
        D = a55;
    }

    return std::pow(eps, m) * D * n * wavenumber *
           std::cosh(n * wavenumber * (waterdepth + (z - zsl))) *
           std::cos(n * phase);
}

// Based on Fenton 1985
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void stokes_waves(
    int StokesOrder,
    amrex::Real wavelength,
    amrex::Real waterdepth,
    amrex::Real waveheight,
    amrex::Real zsl,
    amrex::Real g,
    amrex::Real x,
    amrex::Real z,
    amrex::Real time,
    amrex::Real& eta,
    amrex::Real& u_w,
    amrex::Real& v_w,
    amrex::Real& w_w)
{
    const amrex::Real wavenumber = 2.0 * M_PI / wavelength;

    // some parameters
    amrex::Real c0{0.0};
    amrex::Real a11{0.0}, a22{0.0}, b22{0.0}, c2{0.0}, d2{0.0}, e2{0.0};
    amrex::Real a31{0.0}, a33{0.0}, b31{0.0}, a42{0.0}, a44{0.0};
    amrex::Real b42{0.0}, b44{0.0}, c4{0.0}, d4{0.0}, e4{0.0};
    amrex::Real a51{0.0}, a53{0.0}, a55{0.0}, b53{0.0}, b55{0.0};

    stokes_coefficients(
        StokesOrder, wavenumber, waterdepth, c0, a11, a22, b22, c2, d2, e2, a31,
        a33, b31, a42, a44, b42, b44, c4, d4, e4, a51, a53, a55, b53, b55);

    const amrex::Real eps = wavenumber * waveheight / 2.; // Steepness (ka)
    const amrex::Real c = (c0 + std::pow(eps, 2) * c2 + std::pow(eps, 4) * c4) *
                          std::sqrt(g / wavenumber);
    // const amrex::Real Q =
    //    c * waterdepth * std::sqrt(std::pow(k, 3) / g) +
    //    d2 * std::pow(eps, 2) +
    //    d4 * std::pow(eps, 4) * std::sqrt(g / std::pow(k, 3));

    const amrex::Real omega = c * wavenumber;
    const amrex::Real phase = wavenumber * x - omega * time;

    eta =
        (eps * std::cos(phase)                           // first order term
         + std::pow(eps, 2) * b22 * std::cos(2. * phase) // second order term
         + std::pow(eps, 3) * b31 * (std::cos(phase) - std::cos(3. * phase)) +
         std::pow(eps, 4) * b42 *
             (std::cos(2. * phase) + b44 * std::cos(4. * phase)) +
         std::pow(eps, 5) *
             (-(b53 + b55) * std::cos(phase) + b53 * std::cos(3 * phase) +
              b55 * std::cos(5 * phase))) + zsl;

    // Compute cosines
    const amrex::Real cc11 = stokes_cosh_cos(
        1, 1, phase, a11, a22, a31, a33, a42, a44, a51, a53, a55, eps,
        wavenumber, waterdepth, zsl, z);
    const amrex::Real cc22 = stokes_cosh_cos(
        2, 2, phase, a11, a22, a31, a33, a42, a44, a51, a53, a55, eps,
        wavenumber, waterdepth, zsl, z);
    const amrex::Real cc31 = stokes_cosh_cos(
        3, 1, phase, a11, a22, a31, a33, a42, a44, a51, a53, a55, eps,
        wavenumber, waterdepth, zsl, z);
    const amrex::Real cc33 = stokes_cosh_cos(
        3, 3, phase, a11, a22, a31, a33, a42, a44, a51, a53, a55, eps,
        wavenumber, waterdepth, zsl, z);
    const amrex::Real cc42 = stokes_cosh_cos(
        4, 2, phase, a11, a22, a31, a33, a42, a44, a51, a53, a55, eps,
        wavenumber, waterdepth, zsl, z);
    const amrex::Real cc44 = stokes_cosh_cos(
        4, 4, phase, a11, a22, a31, a33, a42, a44, a51, a53, a55, eps,
        wavenumber, waterdepth, zsl, z);
    const amrex::Real cc51 = stokes_cosh_cos(
        5, 1, phase, a11, a22, a31, a33, a42, a44, a51, a53, a55, eps,
        wavenumber, waterdepth, zsl, z);
    const amrex::Real cc53 = stokes_cosh_cos(
        5, 3, phase, a11, a22, a31, a33, a42, a44, a51, a53, a55, eps,
        wavenumber, waterdepth, zsl, z);
    const amrex::Real cc55 = stokes_cosh_cos(
        5, 5, phase, a11, a22, a31, a33, a42, a44, a51, a53, a55, eps,
        wavenumber, waterdepth, zsl, z);

    // Compute sines
    const amrex::Real ss11 = stokes_sinh_sin(
        1, 1, phase, a11, a22, a31, a33, a42, a44, a51, a53, a55, eps,
        wavenumber, waterdepth, zsl, z);
    const amrex::Real ss22 = stokes_sinh_sin(
        2, 2, phase, a11, a22, a31, a33, a42, a44, a51, a53, a55, eps,
        wavenumber, waterdepth, zsl, z);
    const amrex::Real ss31 = stokes_sinh_sin(
        3, 1, phase, a11, a22, a31, a33, a42, a44, a51, a53, a55, eps,
        wavenumber, waterdepth, zsl, z);
    const amrex::Real ss33 = stokes_sinh_sin(
        3, 3, phase, a11, a22, a31, a33, a42, a44, a51, a53, a55, eps,
        wavenumber, waterdepth, zsl, z);
    const amrex::Real ss42 = stokes_sinh_sin(
        4, 2, phase, a11, a22, a31, a33, a42, a44, a51, a53, a55, eps,
        wavenumber, waterdepth, zsl, z);
    const amrex::Real ss44 = stokes_sinh_sin(
        4, 4, phase, a11, a22, a31, a33, a42, a44, a51, a53, a55, eps,
        wavenumber, waterdepth, zsl, z);
    const amrex::Real ss51 = stokes_sinh_sin(
        5, 1, phase, a11, a22, a31, a33, a42, a44, a51, a53, a55, eps,
        wavenumber, waterdepth, zsl, z);
    const amrex::Real ss53 = stokes_sinh_sin(
        5, 3, phase, a11, a22, a31, a33, a42, a44, a51, a53, a55, eps,
        wavenumber, waterdepth, zsl, z);
    const amrex::Real ss55 = stokes_sinh_sin(
        5, 5, phase, a11, a22, a31, a33, a42, a44, a51, a53, a55, eps,
        wavenumber, waterdepth, zsl, z);

    u_w = c0 * std::sqrt(g / std::pow(wavenumber, 3)) *
          (cc11 + cc22 + cc31 + cc33 + cc42 + cc44 + cc51 + cc53 + cc55);
    v_w = 0.0;
    w_w = c0 * std::sqrt(g / std::pow(wavenumber, 3)) *
          (ss11 + ss22 + ss31 + ss33 + ss42 + ss44 + ss51 + ss53 + ss55);
}

} // namespace amr_wind::ocean_waves::relaxation_zones
#endif
