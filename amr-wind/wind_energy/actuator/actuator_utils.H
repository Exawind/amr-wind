#ifndef ACUTATOR_UTILS_H
#define ACUTATOR_UTILS_H

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
 *  \param rbox The bounding box that defines the region of influence of a turbine
 */
std::set<int> determine_influenced_procs(
    const amrex::AmrCore& mesh, const amrex::RealBox& rbox);

} // namespace utils
} // namespace actuator
} // namespace amr_wind

#endif /* ACUTATOR_UTILS_H */
