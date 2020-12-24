#ifndef WING_OPS_H
#define WING_OPS_H

#include "amr-wind/wind_energy/actuator/wing/ActuatorWing.H"
#include "amr-wind/wind_energy/actuator/actuator_types.H"

namespace amr_wind {
namespace actuator {
namespace utils {
class ActParser;
}

namespace wing {

/** Parse common data for modeling fixed wings as actuators
 */
void read_inputs(WingBaseData&, ActInfo&, const utils::ActParser&);

/** Initialize core data structures when modeling fixed wings
 */
void init_data_structures(WingBaseData&, ActGrid&);

} // namespace wing
} // namespace actuator
} // namespace amr_wind

#endif /* WING_OPS_H */
