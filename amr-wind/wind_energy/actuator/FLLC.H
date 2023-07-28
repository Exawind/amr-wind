#ifndef FLLC
#define FLLC

#include "amr-wind/wind_energy/actuator/actuator_types.H"
#include "amr-wind/wind_energy/actuator/ActParser.H"

namespace amr_wind::actuator {

enum class FLLCType { ConstantChord, VariableChord };

static std::map<std::string, FLLCType> FLLCTypeMap{
    {"constant_chord", FLLCType::ConstantChord},
    {"variable_chord", FLLCType::VariableChord}};

struct FLLCData
{
    // constants
    FLLCType correction_type{FLLCType::VariableChord};
    amrex::Real epsilon;
    amrex::Real relaxation_factor{0.1};
    amrex::Real fllc_start_time{0.0};
    RealList dr;
    RealList optimal_epsilon;

    // computed values
    VecList force_point_velocity;
    VecList les_velocity;
    VecList optimal_velocity;
    VecList correction_velocity;
    VecList lift;
    VecList grad_lift;
    bool initialized{false};
    bool different_sizes;
    RealList span_distance_vel;
    RealList span_distance_force;
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

} // namespace amr_wind::actuator

#endif /* FLLC */
