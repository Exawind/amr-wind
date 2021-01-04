#ifndef TURBINEFAST_H
#define TURBINEFAST_H

#include "amr-wind/wind_energy/actuator/turbine/turbine_types.H"
#include "amr-wind/wind_energy/actuator/turbine/fast/fast_types.H"
#include "amr-wind/wind_energy/actuator/turbine/fast/FastIface.H"
#include "amr-wind/core/ExtSolver.H"

namespace amr_wind {
namespace actuator {

struct TurbineFastData : public TurbineBaseData
{
    amrex::Real density{1.0};

    ::exw_fast::FastTurbine fast_data;
    ::exw_fast::FastIface* fast{nullptr};

    MPI_Comm tcomm{MPI_COMM_NULL};
};

struct TurbineFast : public TurbineType
{
    using InfoType = TurbineInfo;
    using GridType = ActGrid;
    using MetaType = TurbineFastData;
    using DataType = ActDataHolder<TurbineFast>;

    static const std::string identifier() { return "TurbineFast"; }
};

} // namespace actuator
} // namespace amr_wind

#endif /* TURBINEFAST_H */
