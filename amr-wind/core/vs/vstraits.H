#ifndef VSTRAITS_H
#define VSTRAITS_H

#include <limits>
#include "AMReX_REAL.H"

namespace amr_wind::vs {

template <typename T>
struct DTraits
{};

template <>
struct DTraits<int>
{
    static constexpr int zero() noexcept { return 0; }
    static constexpr int one() noexcept { return 1; }
    static constexpr int max() noexcept
    {
        return std::numeric_limits<int>::max();
    }
    static constexpr int min() noexcept
    {
        return std::numeric_limits<int>::min();
    }
};

template <>
struct DTraits<double>
{
    static constexpr double zero() noexcept { return 0.0; }
    static constexpr double one() noexcept { return 1.0; }
    static constexpr double max() noexcept
    {
        return std::numeric_limits<double>::max();
    }
    static constexpr double min() noexcept
    {
        return std::numeric_limits<double>::min();
    }
    static constexpr double eps() noexcept { return 1.0e-15; }
};

template <>
struct DTraits<float>
{
    static constexpr float zero() noexcept { return 0.0F; }
    static constexpr float one() noexcept { return 1.0F; }
    static constexpr float max() noexcept
    {
        return std::numeric_limits<float>::max();
    }
    static constexpr float min() noexcept
    {
        return std::numeric_limits<float>::min();
    }
    static constexpr float eps() noexcept { return 1.0e-7F; }
};

} // namespace amr_wind::vs

#endif /* VSTRAITS_H */
