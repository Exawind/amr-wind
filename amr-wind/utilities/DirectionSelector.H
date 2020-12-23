//
//  DirectionSelector.H
//  amr-wind
//
//

#ifndef DirectionSelector_H
#define DirectionSelector_H

#include "AMReX_Box.H"
#include "AMReX_Gpu.H"

/** select index based on direction input */
template <int Index>
struct DirectionSelector
{
    int operator()(int i, int j, int k) const;
};

template <>
struct DirectionSelector<0>
{
    AMREX_GPU_HOST_DEVICE int operator()(int i, int, int) const { return i; };
};

template <>
struct DirectionSelector<1>
{
    AMREX_GPU_HOST_DEVICE int operator()(int, int j, int) const { return j; };
};

template <>
struct DirectionSelector<2>
{
    AMREX_GPU_HOST_DEVICE int operator()(int, int, int k) const { return k; };
};

using XDir = DirectionSelector<0>;
using YDir = DirectionSelector<1>;
using ZDir = DirectionSelector<2>;

// Given a box, return a 2D box perpendicular to the selected axis.
// For example, if we're using ZDir, return a box covering the x-y plane.
// The IntVect is used to set the constant index in the parallel direction.
template <typename IndexSelector>
AMREX_GPU_HOST_DEVICE amrex::Box
PerpendicularBox(const amrex::Box& bx, const amrex::IntVect& iv)
{
    amrex::IntVect plane_lo, plane_hi;

    if (std::is_same<IndexSelector, XDir>::value) {
        plane_lo = {iv[0], bx.smallEnd(1), bx.smallEnd(2)};
        plane_hi = {iv[0], bx.bigEnd(1), bx.bigEnd(2)};
    } else if (std::is_same<IndexSelector, YDir>::value) {
        plane_lo = {bx.smallEnd(0), iv[1], bx.smallEnd(2)};
        plane_hi = {bx.bigEnd(0), iv[1], bx.bigEnd(2)};
    } else {
        plane_lo = {bx.smallEnd(0), bx.smallEnd(1), iv[2]};
        plane_hi = {bx.bigEnd(0), bx.bigEnd(1), iv[2]};
    }

    amrex::Box pbx(plane_lo, plane_hi);

    return pbx;
}

// Given a box, return a 1D box parallel to the selected axis.
// For example, if we're using ZDir, return a box covering the z axis.
// The IntVect is used to set the constant indices in the perpendicular
// direction.
template <typename IndexSelector>
AMREX_GPU_HOST_DEVICE amrex::Box
ParallelBox(const amrex::Box& bx, const amrex::IntVect& iv)
{
    amrex::IntVect line_lo, line_hi;

    if (std::is_same<IndexSelector, XDir>::value) {
        line_lo = {bx.smallEnd(0), iv[1], iv[2]};
        line_hi = {bx.bigEnd(0), iv[1], iv[2]};
    } else if (std::is_same<IndexSelector, YDir>::value) {
        line_lo = {iv[0], bx.smallEnd(1), iv[2]};
        line_hi = {iv[0], bx.bigEnd(1), iv[2]};
    } else {
        line_lo = {iv[0], iv[1], bx.smallEnd(2)};
        line_hi = {iv[0], iv[1], bx.bigEnd(2)};
    }

    amrex::Box lbx(line_lo, line_hi);

    return lbx;
}

#endif /* DirectionSelector_H */
