#ifndef FLAT_H
#define FLAT_H

#include "amr-wind/immersed_boundary/complex_terrain/ComplexTerrain.H"

namespace amr_wind {
namespace ib {

struct FlatData : public ComplexTerrainBaseData
{
    amrex::Real height{0.04};
};

struct Flat : public ComplexTerrainType
{
    using InfoType = IBInfo;
    using MetaType = FlatData;
    using DataType = IBDataHolder<Flat>;

    static const std::string identifier() { return "Flat"; }
};

} // namespace ib
} // namespace amr_wind

#endif /* Flat_H */
