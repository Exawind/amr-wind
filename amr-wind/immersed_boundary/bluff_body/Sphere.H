#ifndef SPHERE_H
#define SPHERE_H

#include "amr-wind/immersed_boundary/bluff_body/BluffBody.H"

namespace amr_wind::ib {

struct SphereData : public BluffBodyBaseData
{
    vs::Vector center_loc{0.0, 0.0, 0.0};
    amrex::Real radius{0.1};
};

struct Sphere : public BluffBodyType
{
    using InfoType = IBInfo;
    using MetaType = SphereData;
    using DataType = IBDataHolder<Sphere>;

    static std::string identifier() { return "Sphere"; }
};

} // namespace amr_wind::ib

#endif /* Sphere_H */
