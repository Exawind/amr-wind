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
    std::unique_ptr<kynema::interfaces::TurbineInterface> interface;

    /** Data access functions needed for every turbine type */

    float* position_at_vel(int dir) const override { return nullptr; }

    float* solid_velocity(int dir) const override { return nullptr; }

    float* fluid_velocity(int dir) const override { return nullptr; }

    float* force(int dir) const override { return nullptr; }

    float* position_at_force(int dir) const override { return nullptr; }

    float* chord_at_force() const override { return nullptr; }

    float* orientation() const override { return nullptr; }

    int length_fluid_velocity(int dir) const override {
        return interface->num_aero_sections();
    }

    int length_force(int dir) const override
    {
        int n_nodes{0};
        for (int n = 0; n < interface->Turbine().blades.size(); ++n) {
            n_nodes += interface->Turbine().blades[n].nodes.size();
        }
        n_nodes += interface->Turbine().tower.nodes.size();
        return n_nodes;
    }

    int length_position_at_force(int dir) const override
    {
        return length_force(dir);
    }

    int length_orientation() const override { return 0; }
};

struct KynemaSolverData
{
    std::string m_output_dir{"kynema_velocity_data"};
};

} // namespace ext_turb

#endif /* KYNEMA_TYPES_H */
