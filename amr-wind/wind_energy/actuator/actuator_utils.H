#ifndef ACUTATOR_UTILS_H
#define ACUTATOR_UTILS_H

#include "amr-wind/core/vs/vector.H"
#include "amr-wind/utilities/trig_ops.H"
#include "AMReX_AmrCore.H"

#include <set>

namespace amr_wind {
namespace actuator {
namespace utils {

/** Return a set of process IDs (MPI ranks) that contain AMR boxes that interact
 *  with a given actuator body.
 *
 *  The region of influence of the actuator body is determined by checking for
 *  intersections with a bounding box ``rbox``.
 *
 *  \param mesh AMReX mesh instance
 *  \param rbox The bounding box that defines the region of influence of a
 * turbine
 */
std::set<int> determine_influenced_procs(
    const amrex::AmrCore& mesh, const amrex::RealBox& rbox);

/** Return the Gaussian smearing factor
 *
 *  \param dist Distance vector of the cell center from the actuator node in
 *  local frame of reference.
 *
 *  \param eps Three-dimensional Gaussian scaling factor
 *
 *  \return Gaussian smearing factor
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
gaussian3d(const vs::Vector& dist, const vs::Vector& eps)
{
    // const amrex::Real fac = 1.0 / std::sqrt(pi() * pi() * pi());
    constexpr amrex::Real fac = 0.17958712212516656;
    const amrex::Real eps_fac = eps.x() * eps.y() * eps.z();
    const vs::Vector rr{dist.x() / eps.x(), dist.y() / eps.y(),
                        dist.z() / eps.z()};

    return (fac / eps_fac) *
           std::exp(-rr.x() * rr.x() - rr.y() * rr.y() - rr.z() * rr.z());
}

} // namespace utils
} // namespace actuator
} // namespace amr_wind

#endif /* ACUTATOR_UTILS_H */
