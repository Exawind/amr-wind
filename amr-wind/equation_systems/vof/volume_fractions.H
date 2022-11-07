#ifndef VOLUME_FRACTIONS_H_
#define VOLUME_FRACTIONS_H_

#include <AMReX_FArrayBox.H>
#include <cmath>

namespace amr_wind::multiphase {

/** Young's finite-difference gradient scheme.
 *
 *  The gradient is computed with a multiplicative factor of -32:
 *  \f[ mm = - 32 * \grad (c) \f]
 *
 *  Known problems: the index (1,1,1), i.e. the central cell in the block, never
 *  occurs: Therefore an isolated droplet will have a normal with all components
 *  to zero.
 */

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void youngs_fd_normal(
    int i,
    int j,
    int k,
    amrex::Array4<amrex::Real const> const& volfrac,
    amrex::Real& mx,
    amrex::Real& my,
    amrex::Real& mz) noexcept
{
    amrex::Real mm1, mm2;

    mm1 = volfrac(i - 1, j - 1, k - 1) + volfrac(i - 1, j - 1, k + 1) +
          volfrac(i - 1, j + 1, k - 1) + volfrac(i - 1, j + 1, k + 1) +
          2.0 * (volfrac(i - 1, j - 1, k) + volfrac(i - 1, j + 1, k) +
                 volfrac(i - 1, j, k - 1) + volfrac(i - 1, j, k + 1)) +
          4.0 * volfrac(i - 1, j, k);
    mm2 = volfrac(i + 1, j - 1, k - 1) + volfrac(i + 1, j - 1, k + 1) +
          volfrac(i + 1, j + 1, k - 1) + volfrac(i + 1, j + 1, k + 1) +
          2.0 * (volfrac(i + 1, j - 1, k) + volfrac(i + 1, j + 1, k) +
                 volfrac(i + 1, j, k - 1) + volfrac(i + 1, j, k + 1)) +
          4.0 * volfrac(i + 1, j, k);
    mx = mm1 - mm2;

    mm1 = volfrac(i - 1, j - 1, k - 1) + volfrac(i - 1, j - 1, k + 1) +
          volfrac(i + 1, j - 1, k - 1) + volfrac(i + 1, j - 1, k + 1) +
          2.0 * (volfrac(i - 1, j - 1, k) + volfrac(i + 1, j - 1, k) +
                 volfrac(i, j - 1, k - 1) + volfrac(i, j - 1, k + 1)) +
          4.0 * volfrac(i, j - 1, k);
    mm2 = volfrac(i - 1, j + 1, k - 1) + volfrac(i - 1, j + 1, k + 1) +
          volfrac(i + 1, j + 1, k - 1) + volfrac(i + 1, j + 1, k + 1) +
          2.0 * (volfrac(i - 1, j + 1, k) + volfrac(i + 1, j + 1, k) +
                 volfrac(i, j + 1, k - 1) + volfrac(i, j + 1, k + 1)) +
          4.0 * volfrac(i, j + 1, k);
    my = mm1 - mm2;

    mm1 = volfrac(i - 1, j - 1, k - 1) + volfrac(i - 1, j + 1, k - 1) +
          volfrac(i + 1, j - 1, k - 1) + volfrac(i + 1, j + 1, k - 1) +
          2.0 * (volfrac(i - 1, j, k - 1) + volfrac(i + 1, j, k - 1) +
                 volfrac(i, j - 1, k - 1) + volfrac(i, j + 1, k - 1)) +
          4.0 * volfrac(i, j, k - 1);
    mm2 = volfrac(i - 1, j - 1, k + 1) + volfrac(i - 1, j + 1, k + 1) +
          volfrac(i + 1, j - 1, k + 1) + volfrac(i + 1, j + 1, k + 1) +
          2.0 * (volfrac(i - 1, j, k + 1) + volfrac(i + 1, j, k + 1) +
                 volfrac(i, j - 1, k + 1) + volfrac(i, j + 1, k + 1)) +
          4.0 * volfrac(i, j, k + 1);
    mz = mm1 - mm2;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void mixed_youngs_central_normal(
    int i,
    int j,
    int k,
    amrex::Array4<amrex::Real const> const& volfrac,
    amrex::Real& mx,
    amrex::Real& my,
    amrex::Real& mz) noexcept
{
    amrex::Real m[AMREX_SPACEDIM + 1][AMREX_SPACEDIM];
    amrex::Real m1, m2;
    // write the plane as: sgn(mx) X =  my Y +  mz Z + alpha
    //                           m00 X = m01 Y + m02 Z + alpha
    m1 = volfrac(i - 1, j, k - 1) + volfrac(i - 1, j, k + 1) +
         volfrac(i - 1, j - 1, k) + volfrac(i - 1, j + 1, k) +
         volfrac(i - 1, j, k);
    m2 = volfrac(i + 1, j, k - 1) + volfrac(i + 1, j, k + 1) +
         volfrac(i + 1, j - 1, k) + volfrac(i + 1, j + 1, k) +
         volfrac(i + 1, j, k);

    m[0][0] = (m1 > m2) ? 1.0 : -1.0;

    m1 = volfrac(i - 1, j - 1, k) + volfrac(i + 1, j - 1, k) +
         volfrac(i, j - 1, k);
    m2 = volfrac(i - 1, j + 1, k) + volfrac(i + 1, j + 1, k) +
         volfrac(i, j + 1, k);
    m[0][1] = 0.5 * (m1 - m2);

    m1 = volfrac(i - 1, j, k - 1) + volfrac(i + 1, j, k - 1) +
         volfrac(i, j, k - 1);
    m2 = volfrac(i - 1, j, k + 1) + volfrac(i + 1, j, k + 1) +
         volfrac(i, j, k + 1);
    m[0][2] = 0.5 * (m1 - m2);

    // write the plane as: sgn(my) Y =  mx X +  mz Z + alpha,
    //                          m11 Y = m10 X + m12 Z + alpha.

    m1 = volfrac(i - 1, j - 1, k) + volfrac(i - 1, j + 1, k) +
         volfrac(i - 1, j, k);
    m2 = volfrac(i + 1, j - 1, k) + volfrac(i + 1, j + 1, k) +
         volfrac(i + 1, j, k);
    m[1][0] = 0.5 * (m1 - m2);

    m1 = volfrac(i, j - 1, k - 1) + volfrac(i, j - 1, k + 1) +
         volfrac(i + 1, j - 1, k) + volfrac(i - 1, j - 1, k) +
         volfrac(i, j - 1, k);
    m2 = volfrac(i, j + 1, k - 1) + volfrac(i, j + 1, k + 1) +
         volfrac(i + 1, j + 1, k) + volfrac(i - 1, j + 1, k) +
         volfrac(i, j + 1, k);

    m[1][1] = (m1 > m2) ? 1.0 : -1.0;

    m1 = volfrac(i, j - 1, k - 1) + volfrac(i, j, k - 1) +
         volfrac(i, j + 1, k - 1);
    m2 = volfrac(i, j - 1, k + 1) + volfrac(i, j, k + 1) +
         volfrac(i, j + 1, k + 1);
    m[1][2] = 0.5 * (m1 - m2);

    // write the plane as: sgn(mz) Z =  mx X +  my Y + alpha
    //                          m22 Z = m20 X + m21 Y + alpha

    m1 = volfrac(i - 1, j, k - 1) + volfrac(i - 1, j, k + 1) +
         volfrac(i - 1, j, k);
    m2 = volfrac(i + 1, j, k - 1) + volfrac(i + 1, j, k + 1) +
         volfrac(i + 1, j, k);
    m[2][0] = 0.5 * (m1 - m2);

    m1 = volfrac(i, j - 1, k - 1) + volfrac(i, j - 1, k + 1) +
         volfrac(i, j - 1, k);
    m2 = volfrac(i, j + 1, k - 1) + volfrac(i, j + 1, k + 1) +
         volfrac(i, j + 1, k);
    m[2][1] = 0.5 * (m1 - m2);

    m1 = volfrac(i - 1, j, k - 1) + volfrac(i + 1, j, k - 1) +
         volfrac(i, j - 1, k - 1) + volfrac(i, j + 1, k - 1) +
         volfrac(i, j, k - 1);
    m2 = volfrac(i - 1, j, k + 1) + volfrac(i + 1, j, k + 1) +
         volfrac(i, j - 1, k + 1) + volfrac(i, j + 1, k + 1) +
         volfrac(i, j, k + 1);

    m[2][2] = (m1 > m2) ? 1.0 : -1.0;

    // normalize each set (mx,my,mz): |mx|+|my|+|mz| = 1
    amrex::Real t0, t1, t2;

    t0 = amrex::Math::abs(m[0][0]) + amrex::Math::abs(m[0][1]) +
         amrex::Math::abs(m[0][2]);
    m[0][0] = m[0][0] / t0;
    m[0][1] = m[0][1] / t0;
    m[0][2] = m[0][2] / t0;

    t0 = amrex::Math::abs(m[1][0]) + amrex::Math::abs(m[1][1]) +
         amrex::Math::abs(m[1][2]);
    m[1][0] = m[1][0] / t0;
    m[1][1] = m[1][1] / t0;
    m[1][2] = m[1][2] / t0;

    t0 = amrex::Math::abs(m[2][0]) + amrex::Math::abs(m[2][1]) +
         amrex::Math::abs(m[2][2]);
    m[2][0] = m[2][0] / t0;
    m[2][1] = m[2][1] / t0;
    m[2][2] = m[2][2] / t0;

    // choose among the three central schemes */
    t0 = amrex::Math::abs(m[0][0]);
    t1 = amrex::Math::abs(m[1][1]);
    t2 = amrex::Math::abs(m[2][2]);

    int cn = 0;
    if (t1 > t0) {
        t0 = t1;
        cn = 1;
    }

    if (t2 > t0) {
        cn = 2;
    }

    // Youngs-CIAM scheme */
    youngs_fd_normal(i, j, k, volfrac, m[3][0], m[3][1], m[3][2]);
    // normalize the set (mx,my,mz): |mx|+|my|+|mz| = 1
    constexpr amrex::Real tiny = 1e-20;
    t0 = amrex::Math::abs(m[3][0]) + amrex::Math::abs(m[3][1]) +
         amrex::Math::abs(m[3][2]) + tiny;
    m[3][0] = m[3][0] / t0;
    m[3][1] = m[3][1] / t0;
    m[3][2] = m[3][2] / t0;

    // choose between the previous choice and Youngs-CIAM
    t0 = amrex::Math::abs(m[3][0]);
    t1 = amrex::Math::abs(m[3][1]);
    t2 = amrex::Math::abs(m[3][2]);
    if (t1 > t0) {
        t0 = t1;
    }
    if (t2 > t0) {
        t0 = t2;
    }

    if (amrex::Math::abs(m[cn][cn]) > t0 && t0 > 0.0) {
        // second t0 condition is to ensure nonzero normal magnitude
        cn = 3;
    }

    // components of the normal vector */
    mx = m[cn][0];
    my = m[cn][1];
    mz = m[cn][2];
}

/* Computes alpha: m1*x1 + m2* x2 + m3*x3 = alpha
 * given that m1+m2+m3=1 (m1,m2,m3>0) and the volumetric fraction volF
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real volume_intercept(
    amrex::Real b1, amrex::Real b2, amrex::Real b3, amrex::Real volF) noexcept
{
    using namespace amrex;

    amrex::Real const_tiny = std::numeric_limits<amrex::Real>::epsilon();
    // (1) order coefficients: m1<m2<m3;
    // (2) get ranges: V1<V2<v3;
    // (3) limit ch (0.d0 < ch < 0.5d0);
    // (4) calculate alpha

    amrex::Real m1 = amrex::min(b1, b2);
    amrex::Real m3 = amrex::max(b1, b2);
    amrex::Real m2 = b3;

    amrex::Real tmp;
    if (m2 < m1) {
        tmp = m1;
        m1 = m2;
        m2 = tmp;
    } else if (m2 > m3) {
        tmp = m3;
        m3 = m2;
        m2 = tmp;
    }

    amrex::Real m12 = m1 + m2;
    amrex::Real pr = amrex::max<amrex::Real>(6. * m1 * m2 * m3, const_tiny);
    amrex::Real V1 = m1 * m1 * m1 / pr;
    amrex::Real V2 = V1 + 0.5 * (m2 - m1) / m3;

    amrex::Real mm, V3;
    if (m3 < m12) {
        mm = m3;
        V3 = (m3 * m3 * (3.0 * m12 - m3) + m1 * m1 * (m1 - 3.0 * m3) +
              m2 * m2 * (m2 - 3.0 * m3)) /
             pr;
    } else {
        mm = m12;
        V3 = 0.5 * mm / m3;
    }

    amrex::Real ch = amrex::min<amrex::Real>(volF, 1.0 - volF);

    amrex::Real alpha, p, q, p12, teta, cs;
    if (ch < V1) {
        alpha = std::cbrt(pr * ch); // case (1)
    } else if (ch < V2) {
        alpha = 0.50 * (m1 + std::sqrt(m1 * m1 + 8.0 * m2 * m3 * (ch - V1)));
    } else if (ch < V3) {
        p = 2.0 * m1 * m2;
        q = 1.5 * m1 * m2 * (m12 - 2.0 * m3 * ch);
        p12 = std::sqrt(p);
        teta = std::acos(q / (p * p12 + 1e-20)) / 3.0;
        cs = std::cos(teta);
        alpha = p12 * (std::sqrt(3.0 * (1.0 - cs * cs)) - cs) + m12;
    } else if (m12 < m3) {
        alpha = m3 * ch + 0.5 * mm;
    } else {
        p = m1 * (m2 + m3) + m2 * m3 - 0.25;
        q = 1.5 * m1 * m2 * m3 * (0.5 - ch);
        p12 = std::sqrt(p);
        teta = std::acos(q / (p * p12 + 1e-20)) / 3.0;
        cs = std::cos(teta);
        alpha = p12 * (std::sqrt(3.0 * (1.0 - cs * cs)) - cs) + 0.5;
    }

    if (volF > 0.5) {
        alpha = 1.0 - alpha;
    }

    return alpha;
}

/** Computes the "CUT VOLUME" V0 given r0, dr0 and
 *  m_1*x_1 + m_2*x_2 + m_3*x_3 = alpha
 *  (1) move origin to r0 along r ;
 * (2) reflect parallelepiped;
 * (3) limit alpha (0<= al0 <=0.5);
 * (4) order coefficients: b1<b2<b3;
 * (5) calculate volume (NOTE: it is assumed:s0=t0=0; ds0=dt0=1.)
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real cut_volume(
    amrex::Real m1,
    amrex::Real m2,
    amrex::Real m3,
    amrex::Real alpha,
    amrex::Real r0,
    amrex::Real dr0) noexcept
{

    amrex::Real const_tiny = std::numeric_limits<amrex::Real>::epsilon();
    amrex::Real al;

    // move origin to x0
    al = alpha - m1 * r0;

    // reflect the figure when negative coefficients
    al = al + amrex::max<amrex::Real>(0.0, -m1 * dr0) +
         amrex::max<amrex::Real>(0.0, -m2) + amrex::max<amrex::Real>(0.0, -m3);

    // normalized equation: m1*y1 + m2*y2 + m3*y3 = alh, with 0 <= m1 <= m2 <=
    // m3 the problem is then solved again in the unit cube
    amrex::Real tmp = std::abs(m1) * dr0 + std::abs(m2) + std::abs(m3);
    amrex::Real n1 = std::abs(m1) / tmp; // need positive coefficients
    amrex::Real n2 = std::abs(m2) / tmp;
    amrex::Real n3 = std::abs(m3) / tmp;
    al = amrex::max<amrex::Real>(
        0.0, amrex::min<amrex::Real>(
                 1.0, al / tmp)); // get new al within safe limits
    amrex::Real al0 =
        amrex::min<amrex::Real>(al, 1.0 - al); // limit to: 0 < alh < 1/2

    // Order coefficients
    amrex::Real b1 = amrex::min(n1 * dr0, n2); // order coefficients
    amrex::Real b3 = amrex::max(n1 * dr0, n2);
    amrex::Real b2 = n3;

    if (b2 < b1) {
        tmp = b1;
        b1 = b2;
        b2 = tmp;
    } else if (b2 > b3) {
        tmp = b3;
        b3 = b2;
        b2 = tmp;
    }

    amrex::Real b12 = b1 + b2;
    // amrex::Real bm = amrex::min(b12, b3);
    // amrex::Real pr = amrex::max(6.0 * m1 * m2 * m3, const_tiny);

    // Compute volume fraction using PLIC from Scardovelli & Zaleski (JCP 2000),
    // adapted from code in paper by Akio Kawano (Computer & Fluids 2016)
    amrex::Real vm1 = b1;
    amrex::Real vm3 = b3;
    amrex::Real vm2 = b2;
    amrex::Real vm12 = b12;
    amrex::Real a = al0;
    amrex::Real v = 0.0;

    if (a > 0.0) {
        if (a < vm1) {
            v = a * a * a / (6.0 * vm1 * vm2 * vm3);
        } else if (a < vm2) {
            v = a * (a - vm1) / (2.0 * vm2 * vm3) +
                vm1 * vm1 / (6.0 * vm2 * vm3 + const_tiny);
        } else if (a < amrex::min(vm12, vm3)) {
            v = (a * a * (3.0 * vm12 - a) + vm1 * vm1 * (vm1 - 3.0 * a) +
                 vm2 * vm2 * (vm2 - 3.0 * a)) /
                (6.0 * vm1 * vm2 * vm3);
        } else if (vm3 < vm12) {
            v = (a * a * (3.0 - 2.0 * a) + vm1 * vm1 * (vm1 - 3.0 * a) +
                 vm2 * vm2 * (vm2 - 3.0 * a) + vm3 * vm3 * (vm3 - 3.0 * a)) /
                (6.0 * vm1 * vm2 * vm3);
        } else {
            v = (a - 0.5 * vm12) / vm3;
        }
    }

    tmp = v;
    amrex::Real FL3D;
    if (al <= 0.5) {
        FL3D = tmp;
    } else {
        FL3D = (1.0 - tmp);
    }

    return FL3D;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void fit_plane(
    const int i,
    const int j,
    const int k,
    amrex::Array4<amrex::Real const> const& volfrac,
    amrex::Real& mx,
    amrex::Real& my,
    amrex::Real& mz,
    amrex::Real& alpha)
{

    mixed_youngs_central_normal(i, j, k, volfrac, mx, my, mz);

    amrex::Real invx = 1.0;
    amrex::Real invy = 1.0;
    amrex::Real invz = 1.0;

    if (mx < 0.0) {
        mx = -mx;
        invx = -1.0;
    }
    if (my < 0.0) {
        my = -my;
        invy = -1.0;
    }
    if (mz < 0.0) {
        mz = -mz;
        invz = -1.0;
    }

    amrex::Real mm2 = mx + my + mz;
    mx = mx / mm2;
    my = my / mm2;
    mz = mz / mm2;

    alpha = volume_intercept(mx, my, mz, volfrac(i, j, k));

    // Back to the original plane
    mx = invx * mx;
    my = invy * my;
    mz = invz * mz;
    alpha = alpha + amrex::min<amrex::Real>(0.0, mx) +
            amrex::min<amrex::Real>(0.0, my) + amrex::min<amrex::Real>(0.0, mz);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE bool interface_band(
    int i,
    int j,
    int k,
    amrex::Array4<amrex::Real const> const& volfrac) noexcept
{
    constexpr amrex::Real tiny = 1e-12;
    // Check if near interface
    amrex::Real VOF_max = 0.0;
    amrex::Real VOF_min = 1.0;
    bool VOF_mid = false;
    constexpr int n_band = 1; // must be 1, 2, or 3
    for (int ii = -n_band; ii <= n_band; ++ii) {
        for (int jj = -n_band; jj <= n_band; ++jj) {
            for (int kk = -n_band; kk <= n_band; ++kk) {
                amrex::Real VOF = volfrac(i + ii, j + jj, k + kk);
                VOF_max = amrex::max(VOF_max, VOF);
                VOF_min = amrex::min(VOF_min, VOF);
                if (VOF - 1.0 < tiny && VOF > tiny) {
                    VOF_mid = true;
                }
            }
        }
    }
    return (VOF_max != VOF_min || VOF_mid);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real levelset_to_vof(
    int i,
    int j,
    int k,
    amrex::Real eps,
    amrex::Array4<amrex::Real const> const& phi) noexcept
{
    amrex::Real mx, my, mz;
    youngs_fd_normal(i, j, k, phi, mx, my, mz);
    mx = std::abs(mx / 32.);
    my = std::abs(my / 32.);
    mz = std::abs(mz / 32.);
    amrex::Real normL1 = (mx + my + mz);
    mx = mx / normL1;
    my = my / normL1;
    mz = mz / normL1;
    // Make sure that alpha is negative far away from the
    // interface
    amrex::Real alpha;
    if (phi(i, j, k) < -eps) {
        alpha = -1.0;
    } else {
        alpha = phi(i, j, k) / normL1;
        alpha = alpha + 0.5;
    }
    amrex::Real result;
    if (alpha >= 1.0) {
        result = 1.0;
    } else if (alpha <= 0.0) {
        result = 0.0;
    } else {
        result = cut_volume(mx, my, mz, alpha, 0.0, 1.0);
    }

    return result;
}

} // namespace amr_wind::multiphase

#endif // VOLUME_FRACTIONS.H
