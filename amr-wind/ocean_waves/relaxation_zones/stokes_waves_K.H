#ifndef STOKES_WAVES_H_
#define STOKES_WAVES_H_

#include <AMReX_FArrayBox.H>
#include <cmath>

// Stokes waves theory adapted from
// Fenton, J., Fifth Order Stokes Theory for Steady Waves
// Journal of Waterway, Port, Coastal and Ocean Engineering, 1985, 111, 216-234

// Updated Table 1 from Fenton 1985 found in
// Fenton, J., Nonlinear Wave Theories
// The Sea Vol. 9 Ocean Engineering Science, 1990
// https://johndfenton.com/Papers/Fenton90b-Nonlinear-wave-theories.pdf

// Relevant results are summarized in
// Kinnas S. A., Notes on fifth order gravity wave theory
// https://www.caee.utexas.edu/prof/kinnas/ce358/oenotes/kinnas_stokes11.pdf

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
    amrex::Real k = (omega * omega) / g;
    // Cannot skip loop
    amrex::Real f = tol + 1.0;
    while (std::abs(f) > tol && iter < iter_max) {
        // Calculate current constants

        // Exponential definition of S = sech(2kd)
        const amrex::Real S =
            2. * std::exp(2. * k * d) / (std::exp(4. * k * d) + 1.);
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
    int stokes_order,
    amrex::Real wavenumber,
    amrex::Real water_depth,
    amrex::Real& c0,
    amrex::Real& a11,
    amrex::Real& a22,
    amrex::Real& b22,
    amrex::Real& c2,
    amrex::Real& a31,
    amrex::Real& a33,
    amrex::Real& b31,
    amrex::Real& a42,
    amrex::Real& a44,
    amrex::Real& b42,
    amrex::Real& b44,
    amrex::Real& c4,
    amrex::Real& a51,
    amrex::Real& a53,
    amrex::Real& a55,
    amrex::Real& b53,
    amrex::Real& b55)
{

    amrex::Real kd = wavenumber * water_depth;
    if (kd > 50. * M_PI) {
        kd = 50 * M_PI;
    }

    // Exponential definition of S = sech(2kd)
    amrex::Real S = 2. * std::exp(2. * kd) / (std::exp(4. * kd) + 1.);
    amrex::Real C = 1.0 - S;
    amrex::Real Sh = std::sinh(kd);
    amrex::Real Th = std::tanh(kd);
    // Exponential definition of coth(kd)
    amrex::Real CTh = (1. + std::exp(-2 * kd)) / (1. - std::exp(-2 * kd));

    c0 = std::sqrt(Th);
    a11 = 1. / std::sinh(kd); // Hyperbolic cosecant
    // Second order coefficients
    a22 = 3. * std::pow(S, 2) / (2 * std::pow(C, 2));
    b22 = CTh * (1 + 2. * S) / (2 * C);
    c2 = std::sqrt(Th) * (2 + 7 * std::pow(S, 2)) / (4 * std::pow(C, 2));
    if (stokes_order == 2) {
        return;
    }

    // Third order coefficients
    a31 = (-4 - 20 * S + 10 * std::pow(S, 2) - 13 * std::pow(S, 3)) /
          (8 * Sh * std::pow(C, 3));
    a33 =
        (-2 * std::pow(S, 2) + 11 * std::pow(S, 3)) / (8 * Sh * std::pow(C, 3));
    b31 = -3 * (1 + 3 * S + 3 * std::pow(S, 2) + 2 * std::pow(S, 3)) /
          (8 * std::pow(C, 3));
    if (stokes_order == 3) {
        return;
    }

    // Fourth order coefficients
    a42 = (12 * S - 14 * std::pow(S, 2) - 264 * std::pow(S, 3) -
           45 * std::pow(S, 4) - 13 * std::pow(S, 5)) /
          (24 * std::pow(C, 5));
    a44 = (10 * std::pow(S, 3) - 174 * std::pow(S, 4) + 291 * std::pow(S, 5) +
           278 * std::pow(S, 6)) /
          (48 * (3 + 2 * S) * std::pow(C, 5));
    b42 = CTh *
          (6 - 26 * S - 182 * std::pow(S, 2) - 204 * std::pow(S, 3) -
           25 * std::pow(S, 4) + 26 * std::pow(S, 5)) /
          (6 * (3 + 2 * S) * std::pow(C, 4));
    b44 = CTh *
          (24 + 92 * S + 122 * std::pow(S, 2) + 66 * std::pow(S, 3) +
           67 * std::pow(S, 4) + 34 * std::pow(S, 5)) /
          (24 * (3 + 2 * S) * std::pow(C, 4));
    c4 = std::sqrt(Th) *
         (4 + 32 * S - 116 * std::pow(S, 2) - 400 * std::pow(S, 3) -
          71 * std::pow(S, 4) + 146 * std::pow(S, 5)) /
         (32 * std::pow(C, 5));
    if (stokes_order == 4) {
        return;
    }

    // Fifth order coefficients
    a51 =
        (-1184 + 32 * S + 13232 * std::pow(S, 2) + 21712 * std::pow(S, 3) +
         20940 * std::pow(S, 4) + 12554 * std::pow(S, 5) -
         500 * std::pow(S, 6) - 3341 * std::pow(S, 7) - 670 * std::pow(S, 8)) /
        (64 * Sh * (3 + 2 * S) * (4 + S) * std::pow(C, 6));
    a53 = (4 * S + 105 * pow(S, 2) + 198 * std::pow(S, 3) -
           1376 * std::pow(S, 4) - 1302 * std::pow(S, 5) -
           117 * std::pow(S, 6) + 58 * std::pow(S, 7)) /
          (32 * Sh * (3 + 2 * S) * std::pow(C, 6));
    a55 =
        (-6 * std::pow(S, 3) + 272 * std::pow(S, 4) - 1552 * std::pow(S, 5) +
         852 * std::pow(S, 6) + 2029 * std::pow(S, 7) + 430 * std::pow(S, 8)) /
        (64 * Sh * (3 + 2 * S) * (4 + S) * std::pow(C, 6));
    b53 = 9 *
          (132 + 17 * S - 2216 * std::pow(S, 2) - 5897 * std::pow(S, 3) -
           6292 * std::pow(S, 4) - 2687 * std::pow(S, 5) +
           194 * std::pow(S, 6) + 467 * std::pow(S, 7) + 82 * std::pow(S, 8)) /
          (128 * (3 + 2 * S) * (4 + S) * std::pow(C, 6));
    b55 =
        5 *
        (300 + 1579 * S + 3176 * std::pow(S, 2) + 2949 * std::pow(S, 3) +
         1188 * std::pow(S, 4) + 675 * std::pow(S, 5) + 1326 * std::pow(S, 6) +
         827 * std::pow(S, 7) + 130 * std::pow(S, 8)) /
        (384 * (3 + 2 * S) * (4 + S) * std::pow(C, 6));
    if (stokes_order == 5) {
        return;
    }

    if (stokes_order > 5 || stokes_order < 2) {
        amrex::Abort(
            "invalid stokes order specified. It should be between 2,3,4 or 5 ");
    }
}

// Based on Fenton 1985
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void stokes_waves(
    int stokes_order,
    amrex::Real wavelength,
    amrex::Real water_depth,
    amrex::Real wave_height,
    amrex::Real zsl,
    amrex::Real g,
    amrex::Real x,
    amrex::Real z,
    amrex::Real time,
    amrex::Real phase_offset,
    amrex::Real& eta,
    amrex::Real& u_w,
    amrex::Real& v_w,
    amrex::Real& w_w)
{
    const amrex::Real wavenumber = 2.0 * M_PI / wavelength;

    // some parameters
    amrex::Real c0{0.0};
    amrex::Real a11{0.0}, a22{0.0}, b22{0.0}, c2{0.0};
    amrex::Real a31{0.0}, a33{0.0}, b31{0.0}, a42{0.0}, a44{0.0};
    amrex::Real b42{0.0}, b44{0.0}, c4{0.0};
    amrex::Real a51{0.0}, a53{0.0}, a55{0.0}, b53{0.0}, b55{0.0};

    stokes_coefficients(
        stokes_order, wavenumber, water_depth, c0, a11, a22, b22, c2, a31, a33,
        b31, a42, a44, b42, b44, c4, a51, a53, a55, b53, b55);

    const amrex::Real eps = wavenumber * wave_height / 2.; // Steepness (ka)
    const amrex::Real c = (c0 + std::pow(eps, 2) * c2 + std::pow(eps, 4) * c4) *
                          std::sqrt(g / wavenumber);

    const amrex::Real omega = c * wavenumber;
    const amrex::Real phase = wavenumber * x - omega * time - phase_offset;

    eta = (eps * std::cos(phase)                           // first order term
           + std::pow(eps, 2) * b22 * std::cos(2. * phase) // second order term
           + std::pow(eps, 3) * b31 * (std::cos(phase) - std::cos(3. * phase)) +
           std::pow(eps, 4) *
               (b42 * std::cos(2. * phase) + b44 * std::cos(4. * phase)) +
           std::pow(eps, 5) *
               (-(b53 + b55) * std::cos(phase) + b53 * std::cos(3 * phase) +
                b55 * std::cos(5 * phase))) /
              wavenumber +
          zsl;

    // Compute velocities components using Eq.(21) Eq.(23) from Kinnas
    // https://www.caee.utexas.edu/prof/kinnas/ce358/oenotes/kinnas_stokes11.pdf
    // Define coefficients using Eq.(19)

    const int MAX_ORDER = 5;
    amrex::GpuArray<amrex::Real, MAX_ORDER> a;
    if (stokes_order == 2) {
        a[0] = a11;
        a[1] = a22;
    }
    if (stokes_order == 3) {
        a[0] = a11 + (eps * eps) * a31;
        a[1] = a22;
        a[2] = a33;
    }
    if (stokes_order == 4) {
        a[0] = a11 + (eps * eps) * a31;
        a[1] = a22 + (eps * eps) * a42;
        a[2] = a33;
        a[3] = a44;
    }
    if (stokes_order == 5) {
        a[0] = a11 + (eps * eps) * a31 + std::pow(eps, 4) * a51;
        a[1] = a22 + (eps * eps) * a42;
        a[2] = a33 + (eps * eps) * a53;
        a[3] = a44;
        a[4] = a55;
    }

    u_w = 0.;
    v_w = 0.;
    w_w = 0.;
    for (int n = 1; n <= stokes_order; ++n) {
        // Upper bound for deep water case
        // This ensure finite values of velocity for large kd's
        amrex::Real nkdz = n * wavenumber * (water_depth + (z - zsl));
        if (nkdz > 50. * M_PI) {
            nkdz = 50. * M_PI;
        }
        u_w += std::pow(eps, n) * n * a[n - 1] * std::cosh(nkdz) *
               std::cos(n * phase);
        w_w += std::pow(eps, n) * n * a[n - 1] * std::sinh(nkdz) *
               std::sin(n * phase);
    }
    u_w *= (c0 * std::sqrt(g / wavenumber));
    w_w *= (c0 * std::sqrt(g / wavenumber));
}

} // namespace amr_wind::ocean_waves::relaxation_zones
#endif
