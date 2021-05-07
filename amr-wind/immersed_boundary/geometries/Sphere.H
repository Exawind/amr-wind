#ifndef GEOMETRY_H
#define GEOMETRY_H

#include "amr-wind/immersed_boundary/geometries/IBGeometry.H

namespace amr_wind {
namespace ib {

struct SphereData : public IBGeometryBaseData
{};

struct Sphere : public IBGeometryType
{
    using InfoType = IBInfo;
    using GridType = IBGrid;
    using MetaType = GeometryData;
    using DataType = IBDataHolder<Geometry>;

    static const std::string identifier() { return "GeometryGenerator"; }
};

} // namespace ib
} // namespace amr_wind

#endif /* GEOMETRY_H */
