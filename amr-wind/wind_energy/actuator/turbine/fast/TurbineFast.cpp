#include "amr-wind/wind_energy/actuator/turbine/fast/TurbineFast.H"
#include "amr-wind/wind_energy/actuator/turbine/fast/turbine_fast_ops.H"
#include "amr-wind/wind_energy/actuator/ActuatorModel.H"

namespace amr_wind::actuator {

template class ActModel<TurbineFast, ActSrcLine>;
template class ActModel<TurbineFast, ActSrcDisk>;

} // namespace amr_wind::actuator
