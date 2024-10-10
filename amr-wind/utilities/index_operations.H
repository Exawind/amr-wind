#ifndef INDEX_OPERATIONS_H
#define INDEX_OPERATIONS_H
#include <AMReX_REAL.H>
#include <AMReX_Vector.H>
#include <AMReX_Box.H>
#include <AMReX_RealBox.H>
#include <AMReX_RealVect.H>
#include <AMReX_Geometry.H>

namespace amr_wind::utils {

//! Return closest index (from lower) of value in vector
AMREX_FORCE_INLINE int
closest_index(const amrex::Vector<amrex::Real>& vec, const amrex::Real value)
{
    auto const it = std::upper_bound(vec.begin(), vec.end(), value);
    AMREX_ALWAYS_ASSERT(it != vec.end());

    const int idx = static_cast<int>(std::distance(vec.begin(), it));
    return std::max(idx - 1, 0);
}

//! Return indices perpendicular to normal
template <typename T = amrex::GpuArray<int, 2>>
AMREX_FORCE_INLINE T perpendicular_idx(const int normal)
{
    switch (normal) {
    case 0:
        return T{1, 2};
    case 1:
        return T{0, 2};
    case 2:
        return T{0, 1};
    default:
        amrex::Abort("Invalid normal value to determine perpendicular indices");
    }
    return T{-1, -1};
}

//! Check if a point is inside a box
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE bool contains(
    const amrex::Box& box,
    const amrex::RealVect& pos,
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& prob_lo,
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& dxinv)
{
    const amrex::IntVect iv(AMREX_D_DECL(
        static_cast<int>(amrex::Math::floor((pos[0] - prob_lo[0]) * dxinv[0])),
        static_cast<int>(amrex::Math::floor((pos[1] - prob_lo[1]) * dxinv[1])),
        static_cast<int>(
            amrex::Math::floor((pos[2] - prob_lo[2]) * dxinv[2]))));
    return box.contains(iv);
}

/** Convert a bounding box into amrex::Box index space at a given level
 *
 *  \param rbx Bounding box as defined in global domain coordinates
 *  \param geom AMReX geometry information for a given level
 *  \return The Box instance that defines the index space equivalent to bounding
 * boxt
 */
amrex::Box
realbox_to_box(const amrex::RealBox& rbx, const amrex::Geometry& geom);

} // namespace amr_wind::utils

#endif /* INDEX_OPERATIONS_H */
