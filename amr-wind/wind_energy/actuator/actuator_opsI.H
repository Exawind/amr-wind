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

    info.root_proc = *std::min_element(info.procs.begin(), info.procs.end());

    const int iproc = amrex::ParallelDescriptor::MyProc();
    auto in_proc = info.procs.find(iproc);
    info.actuator_in_proc = (in_proc != info.procs.end());
}

} // namespace ops
} // namespace actuator
} // namespace amr_wind

#endif /* ACTUATOR_OPSI_H */
