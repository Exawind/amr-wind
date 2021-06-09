#include "amr-wind/wind_energy/actuator/disk/ConstantCt.H"
#include "amr-wind/wind_energy/actuator/disk/constant_ct_ops.H"
#include "amr-wind/wind_energy/actuator/ActuatorModel.H"

namespace amr_wind {
namespace actuator {

template class ActModel<ConstantCt, ActSrcDisk>;

}
} // namespace amr_wind
