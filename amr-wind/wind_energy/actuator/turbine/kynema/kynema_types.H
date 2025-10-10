#ifndef KYNEMA_TYPES_H
#define KYNEMA_TYPES_H

#include "amr-wind/wind_energy/actuator/turbine/kynema/kynema_wrapper.H"
#include "amr-wind/wind_energy/actuator/turbine/external/external_base_types.H"

#include <string>

namespace ext_turb {

/** Representation of a turbine for exchanging data with Kynema
 */
struct KynemaTurbine : public ext_turb::ExternalTurbine
{
    /** Variables exclusive to Kynema turbine type
     */
    kynema::interfaces::TurbineInterface* interface;

    /** Data access functions needed for every turbine type */

    float* position_at_vel(int dir) const override { return nullptr; }

    float* solid_velocity(int dir) const override { return nullptr; }

    float* fluid_velocity(int dir) const override { return nullptr; }

    float* force(int dir) const override { return nullptr; }

    float* position_at_force(int dir) const override { return nullptr; }

    float* chord_at_force() const override { return nullptr; }

    float* orientation() const override { return nullptr; }

    int length_fluid_velocity(int dir) const override { return 0; }

    int length_force(int dir) const override { return 0; }

    int length_position_at_force(int dir) const override { return 0; }

    int length_orientation() const override { return 0; }
};

struct KynemaSolverData
{
    std::string m_output_dir{"kynema_velocity_data"};
};

} // namespace ext_turb

#endif /* KYNEMA_TYPES_H */
