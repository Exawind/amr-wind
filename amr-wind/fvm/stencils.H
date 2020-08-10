#ifndef STENCILS_H
#define STENCILS_H

/** \file stencils.H
 *  \brief Second-order finite volume discretization stencil coefficients
 *
 *  The stencils are defined using 9 coefficients corresponding to `i+1`, `i`
 *  and `i-1` respectively in each of the three directions (see table below).
 *
 *  Coefficients for the first derivative
 *  ```
 *  |   | +1  | 0   | -1  |
 *  |---+-----+-----+-----|
 *  | i | c00 | c01 | c02 |
 *  | j | c10 | c11 | c12 |
 *  | k | c20 | c21 | c22 |
 *  ```
 */

#include "AMReX_REAL.H"
#include "AMReX_Box.H"
#include "AMReX_Orientation.H"
#include "AMReX_Geometry.H"

namespace amr_wind {
namespace fvm {
namespace stencil {

namespace impl {
inline amrex::Box
box_lo(const amrex::Box& bx, const amrex::Geometry& geom, const int idir)
{
    if (geom.isPeriodic(idir)) return amrex::Box();
    if (bx.smallEnd(idir) != geom.Domain().smallEnd(idir)) return amrex::Box();

    amrex::IntVect low(bx.smallEnd());
    amrex::IntVect hi(bx.bigEnd());
    int sm = low[idir];
    low.setVal(idir, sm);
    hi.setVal(idir, sm);
    return amrex::Box(low, hi);
}

inline amrex::Box
box_hi(const amrex::Box& bx, const amrex::Geometry& geom, const int idir)
{
    if (geom.isPeriodic(idir)) return amrex::Box();
    if (bx.bigEnd(idir) != geom.Domain().bigEnd(idir)) return amrex::Box();

    amrex::IntVect low(bx.smallEnd());
    amrex::IntVect hi(bx.bigEnd());
    int sm = hi[idir];
    low.setVal(idir, sm);
    hi.setVal(idir, sm);
    return amrex::Box(low, hi);
}
} // namespace impl

struct StencilInterior
{
    // First derivatives
    static constexpr amrex::Real c00 =  0.5;
    static constexpr amrex::Real c01 =  0.0;
    static constexpr amrex::Real c02 = -0.5;
    static constexpr amrex::Real c10 =  c00;
    static constexpr amrex::Real c11 =  c01;
    static constexpr amrex::Real c12 =  c02;
    static constexpr amrex::Real c20 =  c00;
    static constexpr amrex::Real c21 =  c01;
    static constexpr amrex::Real c22 =  c02;

    // Second derivatives
    static constexpr amrex::Real s00 =  1.0;
    static constexpr amrex::Real s01 = -2.0;
    static constexpr amrex::Real s02 =  1.0;
    static constexpr amrex::Real s10 =  s00;
    static constexpr amrex::Real s11 =  s01;
    static constexpr amrex::Real s12 =  s02;
    static constexpr amrex::Real s20 =  s00;
    static constexpr amrex::Real s21 =  s01;
    static constexpr amrex::Real s22 =  s02;

    static amrex::Box box(const amrex::Box& bx, const amrex::Geometry&)
    {
        return bx;
    }
};

struct StencilILO
{
    static constexpr int direction = 0;
    static constexpr int side = amrex::Orientation::Side::low;

    static constexpr amrex::Real c00 =  0.5;
    static constexpr amrex::Real c01 =  0.5;
    static constexpr amrex::Real c02 = -1.0;
    static constexpr amrex::Real c10 =  0.5;
    static constexpr amrex::Real c11 =  0.0;
    static constexpr amrex::Real c12 = -0.5;
    static constexpr amrex::Real c20 =  c10;
    static constexpr amrex::Real c21 =  c11;
    static constexpr amrex::Real c22 =  c12;

    static constexpr amrex::Real s00 =  4./3.;  // axx
    static constexpr amrex::Real s01 = -4;      // bxx
    static constexpr amrex::Real s02 =  8./3.;  // cxx
    static constexpr amrex::Real s10 =  1;      // ayy
    static constexpr amrex::Real s11 =  -2;     // byy
    static constexpr amrex::Real s12 =  1;      // cyy
    static constexpr amrex::Real s20 =  s10;    // azz
    static constexpr amrex::Real s21 =  s11;    // bzz
    static constexpr amrex::Real s22 =  s12;    // czz

    static amrex::Box box(const amrex::Box& bx, const amrex::Geometry& geom)
    {
        return impl::box_lo(bx, geom, StencilILO::direction);
    }
};

struct StencilJLO
{
    static constexpr int direction = 1;
    static constexpr int side = amrex::Orientation::Side::low;

    static constexpr amrex::Real c00 =  0.5;
    static constexpr amrex::Real c01 =  0.0;
    static constexpr amrex::Real c02 = -0.5;
    static constexpr amrex::Real c10 =  0.5;
    static constexpr amrex::Real c11 =  0.5;
    static constexpr amrex::Real c12 = -1.0;
    static constexpr amrex::Real c20 =  c00;
    static constexpr amrex::Real c21 =  c01;
    static constexpr amrex::Real c22 =  c02;

    static constexpr amrex::Real s00 =  1.0;    // axx
    static constexpr amrex::Real s01 = -2.0;    // bxx
    static constexpr amrex::Real s02 =  1.0;    // cxx
    static constexpr amrex::Real s10 =  4./3.;  // ayy
    static constexpr amrex::Real s11 =  -4.;    // byy
    static constexpr amrex::Real s12 =  8./3.;  // cyy
    static constexpr amrex::Real s20 =  s00;    // azz
    static constexpr amrex::Real s21 =  s01;    // bzz
    static constexpr amrex::Real s22 =  s02;    // czz

    static amrex::Box box(const amrex::Box& bx, const amrex::Geometry& geom)
    {
        return impl::box_lo(bx, geom, StencilJLO::direction);
    }
};

struct StencilKLO
{
    static constexpr int direction = 2;
    static constexpr int side = amrex::Orientation::Side::low;

    static constexpr amrex::Real c00 =  0.5;
    static constexpr amrex::Real c01 =  0.0;
    static constexpr amrex::Real c02 = -0.5;
    static constexpr amrex::Real c10 =  c00;
    static constexpr amrex::Real c11 =  c01;
    static constexpr amrex::Real c12 =  c02;
    static constexpr amrex::Real c20 =  0.5;
    static constexpr amrex::Real c21 =  0.5;
    static constexpr amrex::Real c22 = -1.0;

    static constexpr amrex::Real s00 =  1.0;    // axx
    static constexpr amrex::Real s01 = -2.0;    // bxx
    static constexpr amrex::Real s02 =  1.0;    // cxx
    static constexpr amrex::Real s10 =  s00;    // ayy
    static constexpr amrex::Real s11 =  s01;    // byy
    static constexpr amrex::Real s12 =  s02;    // cyy
    static constexpr amrex::Real s20 =  4./3.;  // azz
    static constexpr amrex::Real s21 =  -4.;    // bzz
    static constexpr amrex::Real s22 =  8./3.;  // czz

    static amrex::Box box(const amrex::Box& bx, const amrex::Geometry& geom)
    {
        return impl::box_lo(bx, geom, StencilKLO::direction);
    }
};

struct StencilIHI
{
    static constexpr int direction = 0;
    static constexpr int side = amrex::Orientation::Side::high;

    static constexpr amrex::Real c00 =  1.0;
    static constexpr amrex::Real c01 = -0.5;
    static constexpr amrex::Real c02 = -0.5;
    static constexpr amrex::Real c10 =  0.5;
    static constexpr amrex::Real c11 =  0.0;
    static constexpr amrex::Real c12 = -0.5;
    static constexpr amrex::Real c20 =  c10;
    static constexpr amrex::Real c21 =  c11;
    static constexpr amrex::Real c22 =  c12;

    static constexpr amrex::Real s00 = 8./3.;   // axx
    static constexpr amrex::Real s01 = -4.;     // bxx
    static constexpr amrex::Real s02 = 4./3.;   // cxx
    static constexpr amrex::Real s10 =  1.;     // ayy
    static constexpr amrex::Real s11 = -2.;     // byy
    static constexpr amrex::Real s12 =  1.;     // cyy
    static constexpr amrex::Real s20 =  s10;    // azz
    static constexpr amrex::Real s21 =  s11;    // bzz
    static constexpr amrex::Real s22 =  s12;    // czz

    static amrex::Box box(const amrex::Box& bx, const amrex::Geometry& geom)
    {
        return impl::box_hi(bx, geom, StencilIHI::direction);
    }
};

struct StencilJHI
{
    static constexpr int direction = 1;
    static constexpr int side = amrex::Orientation::Side::high;

    static constexpr amrex::Real c00 =  0.5;
    static constexpr amrex::Real c01 =  0.0;
    static constexpr amrex::Real c02 = -0.5;
    static constexpr amrex::Real c10 =  1.0;
    static constexpr amrex::Real c11 = -0.5;
    static constexpr amrex::Real c12 = -0.5;
    static constexpr amrex::Real c20 =  c00;
    static constexpr amrex::Real c21 =  c01;
    static constexpr amrex::Real c22 =  c02;

    static constexpr amrex::Real s00 =  1.0;    // axx
    static constexpr amrex::Real s01 = -2.0;    // bxx
    static constexpr amrex::Real s02 =  1.0;    // cxx
    static constexpr amrex::Real s10 =  8./3.;  // ayy
    static constexpr amrex::Real s11 =  -4.;    // byy
    static constexpr amrex::Real s12 =  4./3.;  // cyy
    static constexpr amrex::Real s20 =  s00;    // azz
    static constexpr amrex::Real s21 =  s01;    // bzz
    static constexpr amrex::Real s22 =  s02;    // czz

    static amrex::Box box(const amrex::Box& bx, const amrex::Geometry& geom)
    {
        return impl::box_hi(bx, geom, StencilJHI::direction);
    }
};

struct StencilKHI
{
    static constexpr int direction = 2;
    static constexpr int side = amrex::Orientation::Side::high;

    static constexpr amrex::Real c00 =  0.5;
    static constexpr amrex::Real c01 =  0.0;
    static constexpr amrex::Real c02 = -0.5;
    static constexpr amrex::Real c10 =  c00;
    static constexpr amrex::Real c11 =  c01;
    static constexpr amrex::Real c12 =  c02;
    static constexpr amrex::Real c20 =  1.0;
    static constexpr amrex::Real c21 = -0.5;
    static constexpr amrex::Real c22 = -0.5;

    static constexpr amrex::Real s00 =  1.0;   // axx
    static constexpr amrex::Real s01 = -2.0;   // bxx
    static constexpr amrex::Real s02 =  1.0;   // cxx
    static constexpr amrex::Real s10 =  s00;   // ayy
    static constexpr amrex::Real s11 =  s01;   // byy
    static constexpr amrex::Real s12 =  s02;   // cyy
    static constexpr amrex::Real s20 =  8./3.; // azz
    static constexpr amrex::Real s21 =  -4.;   // bzz
    static constexpr amrex::Real s22 =  4./3.; // czz

    static amrex::Box box(const amrex::Box& bx, const amrex::Geometry& geom)
    {
        return impl::box_hi(bx, geom, StencilKHI::direction);
    }
};

} // namespace stencil
} // namespace fvm
} // namespace amr_wind

#endif /* STENCILS_H */
