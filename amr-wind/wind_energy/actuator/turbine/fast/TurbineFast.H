#ifndef TURBINEFAST_H
#define TURBINEFAST_H

#include "amr-wind/wind_energy/actuator/turbine/external/TurbineExternal.H"
#include "amr-wind/wind_energy/actuator/turbine/fast/fast_types.H"
 
namespace amr_wind::actuator {

using TurbineFastData = TurbineExternalData<::exw_fast::FastTurbine,::exw_fast::FastSolverData>;
using TurbineFast = TurbineFromExtSolver<::exw_fast::FastTurbine,::exw_fast::FastSolverData>;

} // namespace amr_wind::actuator

#endif /* TURBINEFAST_H */
