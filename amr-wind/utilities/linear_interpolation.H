#ifndef LINEAR_INTERPOLATION_H
#define LINEAR_INTERPOLATION_H

#include "AMReX_Gpu.H"

namespace amr_wind::interp {

enum class Limits : int {
    LOWLIM = -2,
    UPLIM = -1,
    VALID = 0,
};

struct Index
{
    int idx;
    Limits lim;
};

template <typename It, typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Index
check_bounds(const It begin, const It end, const T& x)
{
    const int sz = static_cast<int>(end - begin);

    if ((sz < 2) || (x < *begin)) {
        return Index{0, Limits::LOWLIM};
    }
    if (x > *(begin + (sz - 1))) {
        return Index{sz - 1, Limits::UPLIM};
    }

    return Index{0, Limits::VALID};
}

template <typename It, typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Index
bisection_search(const It begin, const It end, const T& x)
{
    auto idx = check_bounds(begin, end, x);
    if ((idx.lim == Limits::LOWLIM) || (idx.lim == Limits::UPLIM)) {
        return idx;
    }

    int il = 0;
    int ir = end - begin;
    const T xl = *begin;

    while ((ir - il) > 1) {
        int mid = (il + ir) >> 1;
        const T xmid = *(begin + mid);

        if ((x - xmid) * (x - xl) <= 0.0) {
            ir = mid;
        } else {
            il = mid;
        }
    }
    idx.idx = il;
    return idx;
}

template <typename It, typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Index
find_index(const It begin, const It end, const T& x, const int hint = 1)
{
    auto idx = check_bounds(begin, end, x);
    if ((idx.lim == Limits::LOWLIM) || (idx.lim == Limits::UPLIM)) {
        return idx;
    }

    for (It it = (begin + hint); it < end; ++it) {
        if (x <= *it) {
            idx.idx = it - begin - 1;
            break;
        }
    }
    return idx;
}

template <typename C1, typename C2>
inline typename std::iterator_traits<C2>::value_type linear(
    const C1 xbegin,
    const C1 xend,
    const C2 yinp,
    const typename std::iterator_traits<C1>::value_type& xout)
{
    using DType1 = typename std::iterator_traits<C1>::value_type;
    const auto idx = bisection_search(xbegin, xend, xout);

    if ((idx.lim == Limits::LOWLIM) || (idx.lim == Limits::UPLIM)) {
        return yinp[idx.idx];
    }

    if (idx.lim == Limits::VALID) {
        static constexpr DType1 eps = 1.0e-8;
        int j = idx.idx;
        const auto denom = (xbegin[j + 1] - xbegin[j]);
        const auto facR = (denom > eps) ? ((xout - xbegin[j]) / denom) : 1.0;
        const auto facL = static_cast<DType1>(1.0) - facR;
        const auto yout = facL * yinp[j] + facR * yinp[j + 1];
        return yout;
    }

    amrex::Abort("idx undefined");
    exit(1);
}

template <typename C1, typename C2>
inline typename C2::value_type
linear(const C1& xinp, const C2& yinp, const typename C1::value_type& xout)
{
    return linear(xinp.data(), (xinp.data() + xinp.size()), yinp.data(), xout);
}

template <typename C1, typename C2>
inline void
linear_monotonic(const C1& xinp, const C2& yinp, const C1& xout, C2& yout)
{
    static constexpr typename C1::value_type eps = 1.0e-8;
    AMREX_ASSERT(xinp.size() == yinp.size());
    AMREX_ASSERT(xout.size() == yout.size());

    int hint = 1;
    int npts = xout.size();
    for (int i = 0; i < npts; ++i) {
        const auto& x = xout[i];
        const auto idx =
            find_index(xinp.data(), (xinp.data() + xinp.size()), x, hint);

        if ((idx.lim == Limits::LOWLIM) || (idx.lim == Limits::UPLIM)) {
            yout[i] = yinp[idx.idx];
        } else if (idx.lim == Limits::VALID) {
            int j = idx.idx;
            const auto denom = (xinp[j + 1] - xinp[j]);
            const auto facR = (denom > eps) ? ((x - xinp[j]) / denom) : 1.0;
            const auto facL = static_cast<typename C1::value_type>(1.0) - facR;
            yout[i] = facL * yinp[j] + facR * yinp[j + 1];
        }
        hint = idx.idx + 1;
    }
}

template <typename C1, typename C2>
inline void linear(const C1& xinp, const C2& yinp, const C1& xout, C2& yout)
{
    AMREX_ASSERT(
        static_cast<amrex::Long>(xinp.size()) ==
        static_cast<amrex::Long>(yinp.size()));
    AMREX_ASSERT(
        static_cast<amrex::Long>(xout.size()) ==
        static_cast<amrex::Long>(yout.size()));

    int npts = xout.size();
    for (int i = 0; i < npts; ++i) {
        yout[i] = linear(xinp, yinp, xout[i]);
    }
}

} // namespace amr_wind::interp

#endif /* LINEAR_INTERPOLATION_H */
