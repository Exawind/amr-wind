#ifndef HOSWAVES_H
#define HOSWAVES_H

#include "amr-wind/ocean_waves/relaxation_zones/RelaxationZones.H"

namespace amr_wind::ocean_waves {

struct HOSWavesData : public RelaxZonesBaseData
{
    // Prefix for HOS files
    std::string HOS_prefix{"HOSGridData"};
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

struct HOSWaves : public RelaxZonesType
{
    using InfoType = OceanWavesInfo;
    using MetaType = HOSWavesData;
    using DataType = OceanWavesDataHolder<HOSWaves>;

    static std::string identifier() { return "HOSWaves"; }
};

} // namespace amr_wind::ocean_waves

#endif /* HOSWAVES_H */
