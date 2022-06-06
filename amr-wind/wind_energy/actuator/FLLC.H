#ifndef FLLC
#define FLLC

#include "amr-wind/wind_energy/actuator/actuator_types.H"
#include "amr-wind/wind_energy/actuator/ActParser.H"

namespace amr_wind {
namespace actuator {
struct FLLCData
{
    // constants
    amrex::Real epsilon;
    amrex::Real relaxation_factor{0.1};
    RealList dx;
    RealList optimal_epsilon;

    // computed values
    VecList les_velocity;
    VecList optimal_velocity;
    VecList correction_velocity;
    VecList lift;
    VecList grad_lift;
    bool initialized{false};
};

/**
 * @brief Function to capture common parsing requirements for the filtered
 * lifting line correction
 *
 * @param pp Parser
 * @param data Data to be populated from parser
 */
void FLLCParse(const utils::ActParser& pp, FLLCData& data);

/**
 * @brief Initialize FLLC data structure. This should be called at the end of
 * the first ComputeForceOp to ensure the data is fully populated
 *
 * @param data FLLC specific data and scratch fields
 * @param view Component view that has references to the grid data
 * @param eps_chord epsilon chord
 */
void FLLCInit(
    FLLCData& data, const ComponentView& view, const amrex::Real eps_chord);

} // namespace actuator
} // namespace amr_wind

#endif /* FLLC */
