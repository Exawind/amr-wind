#ifndef TURBINE_UTILS_H
#define TURBINE_UTILS_H

#include "amr-wind/wind_energy/actuator/turbine/turbine_types.H"
#include "amr-wind/wind_energy/actuator/ActParser.H"

namespace amr_wind {
namespace actuator {
namespace utils {
void read_inputs(
    TurbineBaseData& tdata, TurbineInfo& info, const utils::ActParser& pp);

void prepare_netcdf_file(
    const std::string&,
    const TurbineBaseData&,
    const TurbineInfo&,
    const ActGrid&);

void write_netcdf(
    const std::string&,
    const TurbineBaseData&,
    const TurbineInfo&,
    const ActGrid&,
    const amrex::Real);

} // namespace utils
} // namespace actuator
} // namespace amr_wind

#endif /* TURBINE_UTILS_H */
