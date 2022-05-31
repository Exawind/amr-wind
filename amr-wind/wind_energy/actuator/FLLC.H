#ifndef FLLC
#define FLLC

#include "amr-wind/wind_energy/actuator/actuator_types.H"
#include "amr-wind/wind_energy/actuator/ActParser.H"

namespace amr_wind {
namespace actuator {
struct FLLCData
{
    FLLCData(VecList& ref_vel): relative_velocity(ref_vel){}

    void init_data(const ComponentView& view, const amrex::Real eps_chord);

    // reference fields
    // assuming this is only going to run on host
    // if we ever run on device we should just move this to a straight
    // copy
    VecList& relative_velocity;

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
};

void FLLCParse(const utils::ActParser& pp, FLLCData& data);

} // namespace actuator
} // namespace amr_wind

#endif /* FLLC */
