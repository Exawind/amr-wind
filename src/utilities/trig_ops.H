#ifndef TRIG_OPS_H
#define TRIG_OPS_H

#include <cmath>
#include "AMREX_REAL.H"

namespace amr_wind {
namespace utils {

constexpr amrex::Real pi()
{
    return M_PI;
}

constexpr amrex::Real two_pi()
{
    return 2.0 * M_PI;
}

constexpr amrex::Real half_pi()
{
    return 0.5 * M_PI;
}

amrex::Real radians(const amrex::Real deg_val)
{
    return pi() * deg_val / 180.0;
}

amrex::Real degrees(const amrex::Real rad_val)
{
    return 180.0 * rad_val / pi();
}

}
}

#endif /* TRIG_OPS_H */
