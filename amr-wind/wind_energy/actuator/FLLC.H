#ifndef FLLC
#define FLLC

#include "amr-wind/wind_energy/actuator/actuator_types.H"
#include "amr-wind/wind_energy/actuator/ActParser.H"

namespace amr_wind {
namespace actuator {
struct FLLCData
{
    void init_data(const ComponentView& view, const amrex::Real eps_chord);

    // constants
    amrex::Real epsilon;
    amrex::Real relaxation_factor{0.1};
    RealList dx;
    RealList optimal_epsilon;

    // computed values
    VecList relative_velocity;
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
