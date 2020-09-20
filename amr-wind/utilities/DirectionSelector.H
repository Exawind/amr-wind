//
//  DirectionSelector.H
//  amr-wind
//
//

#ifndef DirectionSelector_H
#define DirectionSelector_H

#include "AMReX_Gpu.H"

/** select index based on direction input */
template<int Index>
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

#endif /* DirectionSelector_H */
