#include "amr-wind/wind_energy/actuator/disk/UniformCt.H"
#include "amr-wind/wind_energy/actuator/disk/uniform_ct_ops.H"
#include "amr-wind/wind_energy/actuator/ActuatorModel.H"

namespace amr_wind {
namespace actuator {

template class ActModel<UniformCt, ActSrcDiskBase<LinearBasis>>;
template class ActModel<UniformCt, ActSrcDiskBase<UniformGaussian>>;

} // namespace actuator
} // namespace amr_wind
