#ifndef ACTUATOR_OPSI_H
#define ACTUATOR_OPSI_H

#include "amr-wind/wind_energy/actuator/actuator_utils.H"
#include "amr-wind/CFDSim.H"

#include <algorithm>

namespace amr_wind {
namespace actuator {
namespace ops {

template <typename T>
void determine_influenced_procs(typename T::DataType& data)
{
    auto& info = data.m_info;
    info.procs =
        utils::determine_influenced_procs(data.m_sim.mesh(), info.bound_box);

    if (info.root_proc > -1) {
        // During regrid, the influenced processes might have changed and might
        // no longer include the root proc. We insert it back to ensure that it
        // is always present on the list.
        info.procs.insert(info.root_proc);
    }

    const int iproc = amrex::ParallelDescriptor::MyProc();
    auto in_proc = info.procs.find(iproc);
    info.actuator_in_proc = (in_proc != info.procs.end());
}

template <typename T>
void determine_root_proc(
    typename T::DataType& data, amrex::Vector<int>& act_proc_count)
{
    auto& info = data.m_info;
    auto& plist = info.procs;
    bool assigned = false;

    plist =
        utils::determine_influenced_procs(data.m_sim.mesh(), info.bound_box);

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
    }
}

} // namespace ops
} // namespace actuator
} // namespace amr_wind

#endif /* ACTUATOR_OPSI_H */
