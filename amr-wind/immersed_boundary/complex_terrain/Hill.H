#ifndef HILL_H
#define HILL_H

#include "amr-wind/immersed_boundary/complex_terrain/ComplexTerrain.H"

namespace amr_wind {
namespace ib {

struct HillData : public ComplexTerrainBaseData
{
    vs::Vector center_loc{0.0, 0.0, 0.0};
    amrex::Real height{0.04};
    amrex::Real half_width{0.1};
};

struct Hill : public ComplexTerrainType
{
    using InfoType = IBInfo;
    using MetaType = HillData;
    using DataType = IBDataHolder<Hill>;

    static const std::string identifier() { return "Hill"; }
};

} // namespace ib
} // namespace amr_wind

#endif /* Hill_H */
