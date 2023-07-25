#ifndef W2AWAVES_H
#define W2AWAVES_H

#include "amr-wind/ocean_waves/relaxation_zones/RelaxationZones.H"
#ifdef AMR_WIND_USE_W2A
#include <Waves2AMR.H>
#endif

namespace amr_wind::ocean_waves {

struct W2AWavesData : public RelaxZonesBaseData
{
    // Prefix for HOS files
    std::string modes_file{"modes_HOS_SWENSE.dat"};
    // Index of first timestep to be used (in filename)
    int HOS_n0{0};
    // Current index
    int HOS_n{0};
    // Time when ow arrays were last changed (last timestep)
    amrex::Real t_last{-1.0};
    // Time of current HOS lists
    amrex::Real HOS_t{0.0};
    // Timestep from HOS data
    amrex::Real HOS_dt{0.0};
};

struct W2AWaves : public RelaxZonesType
{
    using InfoType = OceanWavesInfo;
    using MetaType = W2AWavesData;
    using DataType = OceanWavesDataHolder<W2AWaves>;

    static std::string identifier() { return "W2AWaves"; }
};

} // namespace amr_wind::ocean_waves

#endif /* W2AWAVES_H */
