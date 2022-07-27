#include "amr-wind/wind_energy/actuator/actuator_utils.H"
#include "amr-wind/wind_energy/actuator/actuator_types.H"
#include "amr-wind/core/MeshMap.H"
#include "amr-wind/CFDSim.H"

namespace amr_wind {
namespace actuator {
namespace utils {

namespace {
/** Convert a bounding box into amrex::Box index space at a given level
 *
 *  \param rbx Bounding box as defined in global domain coordinates
 *  \param geom AMReX geometry information for a given level
 *  \param map Pointer to mesh mapping to enable transformation if it is used
 *  \return The Box instance that defines the index space equivalent to bounding
 * boxt
 */
amrex::Box realbox_to_box(
    const amrex::RealBox& rbx,
    const amrex::Geometry& geom,
    const ::amr_wind::MeshMap* map = nullptr)
{
    const auto* problo = geom.ProbLo();
    const auto* probhi = geom.ProbHi();
    const auto* dxi = geom.InvCellSize();
    amrex::Vector<amrex::Real> real_lo, real_hi;
    if (map != nullptr) {
        real_lo = map->unmap(rbx.lo(), geom);
        real_hi = map->unmap(rbx.hi(), geom);
    } else {
        real_lo = {rbx.lo()[0], rbx.lo()[1], rbx.lo()[2]};
        real_hi = {rbx.hi()[0], rbx.hi()[1], rbx.hi()[2]};
    }

    const amrex::RealBox rbx_transformmed(&real_lo[0], &real_hi[0]);

    amrex::IntVect lo, hi;
    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
        amrex::Real bbox_min = amrex::max(rbx_transformmed.lo()[i], problo[i]);
        amrex::Real bbox_max = amrex::min(rbx_transformmed.hi()[i], probhi[i]);

        amrex::Real rlo = amrex::Math::floor((bbox_min - problo[i]) * dxi[i]);
        amrex::Real rhi = amrex::Math::ceil((bbox_max - problo[i]) * dxi[i]);

        lo[i] = static_cast<int>(rlo);
        hi[i] = static_cast<int>(rhi);
    }

    return amrex::Box{lo, hi};
}

} // namespace

std::set<int>
determine_influenced_procs(const CFDSim& sim, const amrex::RealBox& rbx)
{
    std::set<int> procs;
    const auto& mesh = sim.mesh();
    const int finest_level = mesh.finestLevel();
    const int nlevels = mesh.finestLevel() + 1;
    // TODO bx should be mapped to unstretched coordinates here when using mesh
    // mapping, could be implemented inside realbox_to_box as well
    auto bx = realbox_to_box(rbx, mesh.Geom(0), sim.mesh_mapping());

    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& ba = mesh.boxArray(lev);
        const auto& dm = mesh.DistributionMap(lev);

        // Get all possible intersections at this level
        const auto& isects = ba.intersections(bx);

        // Extract the processor ranks
        for (const auto& is : isects) {
            procs.insert(dm[is.first]);
        }

        if (lev < finest_level) {
            bx = bx.refine(mesh.refRatio(lev));
        }
    }

    return procs;
}

void determine_root_proc(ActInfo& info, amrex::Vector<int>& act_proc_count)
{
    auto& plist = info.procs;
    bool assigned = false;

    // If any of the influenced procs is free (i.e., doesn't have a turbine
    // assigned to it) elect it as the root proc for this turbine and return
    // early.
    for (auto ip : plist) {
        if (act_proc_count[ip] < 1) {
            info.root_proc = ip;
            ++act_proc_count[ip];
            assigned = true;
            break;
        }
    }

    // If we found a root proc there is nothing more to do, so return early
    if (assigned) {
        const int iproc = amrex::ParallelDescriptor::MyProc();
        auto in_proc = info.procs.find(iproc);
        info.actuator_in_proc = (in_proc != info.procs.end());
        info.is_root_proc = (info.root_proc == iproc);

        // By default we request all processes where turbine is active to have
        // velocities sampled. Individual actuator instances can override this
        info.sample_vel_in_proc = info.actuator_in_proc;
        return;
    }

    // If we have reached here, then we have more turbines than processes
    // available. We will assign the current turbine to the process that is
    // managing the lowest number of turbines.

    // Determine the MPI rank that contains the fewest turbines
    auto it = std::min_element(act_proc_count.begin(), act_proc_count.end());
    // Make it the root process for this turbine
    info.root_proc = std::distance(act_proc_count.begin(), it);
    // Make sure the root process is part of the process list
    plist.insert(info.root_proc);
    // Increment turbine count with the global tracking array
    ++(*it);

    {
        const int iproc = amrex::ParallelDescriptor::MyProc();
        auto in_proc = info.procs.find(iproc);
        info.actuator_in_proc = (in_proc != info.procs.end());
        info.is_root_proc = (info.root_proc == iproc);

        // By default we request all processes where turbine is active to have
        // velocities sampled. Individual actuator instances can override this
        info.sample_vel_in_proc = info.actuator_in_proc;
    }
}

} // namespace utils
} // namespace actuator
} // namespace amr_wind
