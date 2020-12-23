#include "amr-wind/wind_energy/actuator/actuator_utils.H"

namespace amr_wind {
namespace actuator {
namespace utils {

namespace {

/** Convert a bounding box into amrex::Box index space at a given level
 *
 *  \param rbx Bounding box as defined in global domain coordinates
 *  \param geom AMReX geometry information for a given level
 *  \return The Box instance that defines the index space equivalent to bounding
 * boxt
 */
amrex::Box
realbox_to_box(const amrex::RealBox& rbx, const amrex::Geometry& geom)
{
    const auto* problo = geom.ProbLo();
    const auto* probhi = geom.ProbHi();
    const auto* dxi = geom.InvCellSize();

    amrex::IntVect lo, hi;
    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
        amrex::Real bbox_min = amrex::max(rbx.lo()[i], problo[i]);
        amrex::Real bbox_max = amrex::min(rbx.hi()[i], probhi[i]);

        amrex::Real rlo = amrex::Math::floor((bbox_min - problo[i]) * dxi[i]);
        amrex::Real rhi = amrex::Math::ceil((bbox_max - problo[i]) * dxi[i]);

        lo[i] = static_cast<int>(rlo);
        hi[i] = static_cast<int>(rhi);
    }

    return amrex::Box{lo, hi};
}

} // namespace

std::set<int> determine_influenced_procs(
    const amrex::AmrCore& mesh, const amrex::RealBox& rbx)
{
    std::set<int> procs;
    const int nlevels = mesh.finestLevel() + 1;
    auto bx = realbox_to_box(rbx, mesh.Geom(0));

    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& ba = mesh.boxArray(lev);
        const auto& dm = mesh.DistributionMap(lev);

        // Get all possible intersections at this level
        const auto& isects = ba.intersections(bx);

        // Extract the processor ranks
        for (const auto& is : isects) procs.insert(dm[is.first]);

        bx = bx.refine(mesh.refRatio(lev));
    }

    return procs;
}

} // namespace utils
} // namespace actuator
} // namespace amr_wind
