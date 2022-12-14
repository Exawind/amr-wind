#ifndef ACTUATORWING_H
#define ACTUATORWING_H

#include "amr-wind/wind_energy/actuator/actuator_types.H"
#include "amr-wind/core/vs/vector_space.H"

#include <string>

namespace amr_wind::actuator {

struct FLLCData;

/** Base data representation of a fixed-wing.
 *
 *  \ingroup actuator
 */
struct WingBaseData
{
    //! Number of points along the wing
    int num_pts;

    //! Starting coordinate of the wing
    vs::Vector start;

    //! Ending coordinate of the wing
    vs::Vector end;

    //! Gaussian smearing factor input by user
    vs::Vector eps_inp;

    //! epsilon / chord
    vs::Vector epsilon_chord;

    //! Pitch angle for the wing
    amrex::Real pitch{0.0};

    //! Relative velocity at the actuator node
    VecList vel_rel;

    //! Filtered Lifting Line Correction data
    std::unique_ptr<FLLCData> fllc;

    //! Data structure required for passing data to the FLLC
    ComponentView component_view;

    //! Width of actuator sections
    RealList dx;

    //! Angle of attack at actuator nodes
    RealList aoa;

    //! Computed lift coefficient at actuator nodes
    RealList cl;

    //! Computed drag coefficient at actuator nodes
    RealList cd;

    //! Chord length at the actuator nodes
    RealList chord;

    //! Total integrated lift over the wing
    amrex::Real lift{0.0};

    //! Total integrated drag over the wing
    amrex::Real drag{0.0};
};

struct WingType : public ActuatorType
{};

} // namespace amr_wind::actuator

#endif /* ACTUATORWING_H */
