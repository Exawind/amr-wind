#ifndef ALGOIM_BOUNDINGBOX_HPP
#define AlGOIM_BOUNDINGBOX_HPP

// Algoim::BoundingBox<T,N>

#include "algoim_real.hpp"
#include "algoim_blitzinc.hpp"

namespace Algoim
{
    // BoundingBox<T,N> describes the extent of a hyperrectangle, i.e., the set of points
    // x = (x(0), ..., x(i), ..., x(N-1)) such that xmin(i) <= x(i) <= xmax(i) for all i.
    template<typename T, int N>
    struct BoundingBox
    {
        TinyVector<TinyVector<T,N>,2> range;

        BoundingBox(const TinyVector<T,N>& min, const TinyVector<T,N>& max) : range(min, max) {}

        inline const TinyVector<T,N>& operator() (int i) const
        {
            return range(i);
        }

        inline const TinyVector<T,N>& min() const
        {
            return range(0);
        }

        inline T& min(int i)
        {
            return range(0)(i);
        }

        inline T min(int i) const
        {
            return range(0)(i);
        }

        inline const TinyVector<T,N>& max() const
        {
            return range(1);
        }

        inline T& max(int i)
        {
            return range(1)(i);
        }

        inline T max(int i) const
        {
            return range(1)(i);
        }

        inline TinyVector<T,N> extent() const
        {
            return range(1) - range(0);
        }

        inline T extent(int i) const
        {
            return range(1)(i) - range(0)(i);
        }

        inline TinyVector<Real,N> midpoint() const
        {
            return (range(0) + range(1)) * 0.5;
        }

        inline Real midpoint(int i) const
        {
            return (range(0)(i) + range(1)(i)) * 0.5;
        }

        inline bool operator==(const BoundingBox& x) const
        {
            for (int dim = 0; dim < N; ++dim)
                if (range(0)(dim) != x.range(0)(dim) || range(1)(dim) != x.range(1)(dim))
                    return false;
            return true;
        }
    };
} // namespace Algoim

#endif
