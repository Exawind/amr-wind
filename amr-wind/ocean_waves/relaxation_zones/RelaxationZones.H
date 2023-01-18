#ifndef RELAX_ZONE_H
#define RELAX_ZONE_H

#include "amr-wind/ocean_waves/OceanWavesTypes.H"
#include "amr-wind/core/vs/vector_space.H"

#include <string>

namespace amr_wind::ocean_waves {

/** Base data representation of a regular waves
 *
 *  \ingroup ocean_waves
 */
struct RelaxZonesBaseData
{
    amrex::Real zsl{0.0};

    amrex::Real water_depth{0.5};

    // Wave generation/absorption parameters
    amrex::Real gen_length{4.0};

    amrex::Real beach_length{8.0};

    bool init_wave_field{false};

    bool has_ramp{false};

    bool has_beach{true};
    bool has_outprofile{false};

    amrex::Real ramp_period{2.0};
};

struct RelaxZonesType : public OceanWavesType
{};

} // namespace amr_wind::ocean_waves

#endif // REGULAR_WAVES_H
