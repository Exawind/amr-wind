#ifndef ALGOIM_UTILITY_HPP
#define ALGOIM_UTILITY_HPP

// Minor utility methods used throughout Algoim

#include "algoim_real.hpp"
#include "algoim_blitzinc.hpp"

namespace Algoim
{
    // Square of x
    template<typename T>
    inline T sqr(T x)
    {
        return x*x;
    }

    // Euclidean norm of a vector x
    template<typename T, int N>
    inline T mag(const TinyVector<T,N>& x)
    {
        T res = x(0)*x(0);
        for (int i = 1; i < N; ++i)
            res += x(i)*x(i);
        return std::sqrt(res);
    }

    // Set the dim'th component of x to a given value
    template<typename T, int N>
    inline TinyVector<T,N> setComponent(const TinyVector<T,N>& x, int dim, T value)
    {
        TinyVector<T,N> y = x;
        y(dim) = value;
        return y;
    }

    // Increment the dim'th component of x by a given value
    template<typename T, int N>
    TinyVector<T,N> shift(const TinyVector<T,N>& x, int dim, T value)
    {
        TinyVector<T,N> y = x;
        y(dim) += value;
        return y;
    }
} // namespace Algoim

#endif
