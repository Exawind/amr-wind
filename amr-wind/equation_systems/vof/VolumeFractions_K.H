#ifndef VolumeFraction_K_H_
#define VolumeFraction_K_H_

#include <AMReX_FArrayBox.H>
#include <cmath>

/* Computes alpha: m1*x1 + m2* x2 + m3*x3 = alpha
 * given that m1+m2+m3=1 (m1,m2,m3>0) and the volumetric fraction cc
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
AL3D(amrex::Real mx, amrex::Real my, amrex::Real mz, amrex::Real cc) noexcept
{
    using namespace amrex;

    amrex::Real const_tiny = std::numeric_limits<amrex::Real>::epsilon();
    amrex::Real b1 = mx;
    amrex::Real b2 = my;
    amrex::Real b3 = mz;
    // (1) order coefficients: m1<m2<m3;
    // (2) get ranges: V1<V2<v3;
    // (3) limit ch (0.d0 < ch < 0.5d0);
    // (4) calculate alpha

    amrex::Real m1, m2, m3;
    m1 = std::min(b1, b2);
    m3 = std::max(b1, b2);
    m2 = b3;

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
    amrex::Real pr = std::max(6. * m1 * m2 * m3, const_tiny);
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

    amrex::Real ch = std::min(cc, 1.0 - cc);

    // Original method
    amrex::Real AL3D, p, q, p12, teta, cs;
    if (ch < V1) {
        AL3D = std::pow(pr * ch, 1. / 3.); // case (1)
    } else if (ch < V2) {
        AL3D = 0.50 * (m1 + std::sqrt(m1 * m1 + 8.0 * m2 * m3 * (ch - V1)));
    } else if (ch < V3) {
        p = 2.0 * m1 * m2;
        q = 1.5 * m1 * m2 * (m12 - 2.0 * m3 * ch);
        p12 = std::sqrt(p);
        teta = std::acos(q / (p * p12)) / 3.0;
        cs = std::cos(teta);
        AL3D = p12 * (std::sqrt(3.0 * (1.0 - cs * cs)) - cs) + m12;
    } else if (m12 < m3) {
        AL3D = m3 * ch + 0.5 * mm;
    } else {
        p = m1 * (m2 + m3) + m2 * m3 - 0.25;
        q = 1.5 * m1 * m2 * m3 * (0.5 - ch);
        p12 = std::sqrt(p);
        teta = std::acos(q / (p * p12)) / 3.0;
        cs = std::cos(teta);
        AL3D = p12 * (std::sqrt(3.0 * (1.0 - cs * cs)) - cs) + 0.5;
    }

    if (cc > 0.5) {
        AL3D = 1.0 - AL3D;
    }

    return AL3D;
}

/** Computes the "CUT VOLUME" V0 given r0, dr0 and
 *  m_1*x_1 + m_2*x_2 + m_3*x_3 = alpha
 *  (1) move origin to r0 along r ;
 * (2) reflect parallelepiped;
 * (3) limit alpha (0<= al0 <=0.5);
 * (4) order coefficients: b1<b2<b3;
 * (5) calculate volume (NOTE: it is assumed:s0=t0=0; ds0=dt0=1.)
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real FL3D(
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
    al =
        al + std::max(0.0, -m1 * dr0) + std::max(0.0, -m2) + std::max(0.0, -m3);

    // normalized equation: m1*y1 + m2*y2 + m3*y3 = alh, with 0 <= m1 <= m2 <=
    // m3 the problem is then solved again in the unit cube
    amrex::Real tmp = std::abs(m1) * dr0 + std::abs(m2) + std::abs(m3);
    amrex::Real n1 = std::abs(m1) / tmp; // need positive coefficients
    amrex::Real n2 = std::abs(m2) / tmp;
    amrex::Real n3 = std::abs(m3) / tmp;
    al =
        std::max(0.0, std::min(1.0, al / tmp)); // get new al within safe limits
    amrex::Real al0 = std::min(al, 1.0 - al);   // limit to: 0 < alh < 1/2

    // Order coefficients
    amrex::Real b1, b2, b3;
    b1 = std::min(n1 * dr0, n2); // order coefficients
    b3 = std::max(n1 * dr0, n2);
    b2 = n3;

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
    //amrex::Real bm = std::min(b12, b3);
    //amrex::Real pr = std::max(6.0 * m1 * m2 * m3, const_tiny);

    // Compute volume fraction using Aoki Kawano (Computer & Fluids 2016) method
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
        } else if (a < std::min(vm12, vm3)) {
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
        FL3D = tmp * dr0;
    } else {
        FL3D = (1.0 - tmp) * dr0;
    }

    return FL3D;
}

#endif // VolumeFraction_K.H
