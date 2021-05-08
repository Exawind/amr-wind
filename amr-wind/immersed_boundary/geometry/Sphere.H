#ifndef SPHERE_H
#define SPHERE_H

#include "amr-wind/immersed_boundary/geometry/BluffBody.H"

namespace amr_wind {
namespace ib {

struct SphereData : public BluffBodyBaseData
{
    RealList centre_loc{0.0, 0.0, 0.0};
};

struct Sphere : public BluffBodyType
{
    using InfoType = IBInfo;
    using GridType = IBGrid;
    using MetaType = SphereData;
    using DataType = IBDataHolder<Sphere>;

    static const std::string identifier() { return "Sphere"; }
};

} // namespace ib
} // namespace amr_wind

#endif /* SPHERE_H */
