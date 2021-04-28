#ifndef IBUTILS_H
#define IBUTILS_H

#include "amr-wind/core/vs/vector.H"
#include "amr-wind/utilities/trig_ops.H"
#include "AMReX_AmrCore.H"

#include <set>

namespace amr_wind {
namespace ib {

struct IBInfo;

namespace utils {

/** Return a set of process IDs (MPI ranks) that contain AMR boxes that interact
 *  with a given immersed boundary.
 *
 *  The region of influence of the immersed boundary is determined by checking for
 *  intersections with a bounding box ``rbox``.
 *
 *  \param mesh AMReX mesh instance
 *  \param rbox The bounding box that defines the region of influence of a
 *  turbine
 */
std::set<int> determine_influenced_procs(
    const amrex::AmrCore& mesh, const amrex::RealBox& rbox);

void determine_root_proc(IBInfo&, amrex::Vector<int>&);

/** Evaluate the Dirac delta phi
 *
 *  \param r Normalized distance as defined in https://doi.org/10.1006/jcph.1999.6293
 *
 *  \return Dirac delta phi
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
dirac_delta_phi(const amrex::Real& r)
{
    const amrex::Real abs_r = amrex::Math::abs(r);

    if ((abs_r > 0.5) && (abs_r <= 0.5))
    {
        return 1/6 * (5 - abs_r - std::sqrt(-3 * std::pow((1 - abs_r), 2) + 1));
    }
    else if (abs_r <= 0.5)
    {
        return 1/3 * (1 + std::sqrt(-3 * std::pow(abs_r, 2) + 1));
    }
    else
        return 0.0;
}

/** Evaluate the Dirac delta factor as defined in https://doi.org/10.1006/jcph.1999.6293
 *
 *  \param dist Distance vector of the cell center from the actuator node in
 *  local frame of reference.
 *
 *  \param dx mesh size in each direction
 *
 *  \return Dirac delta factor
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
dirac_delta(const vs::Vector& dist, const vs::Vector& dx)
{
    // clang-format off
    return 1/dx[0] * dirac_delta_phi(dist[0]/dx[0]) +
           1/dx[1] * dirac_delta_phi(dist[1]/dx[1]) +
           1/dx[2] * dirac_delta_phi(dist[2]/dx[2]);
    // clang-format on
}

} // namespace utils
} // namespace ib
} // namespace amr_wind

#endif /* IBUTILS_H */
