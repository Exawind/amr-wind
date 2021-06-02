#ifndef CYLINDER_H
#define CYLINDER_H

#include "amr-wind/immersed_boundary/geometry/BluffBody.H"

namespace amr_wind {
namespace ib {

struct CylinderData : public BluffBodyBaseData
{
    vs::Vector center_loc{0.0, 0.0, 0.0};
    amrex::Real radius{0.1};
    amrex::Real length{1.0};
};

struct Cylinder : public BluffBodyType
{
    using InfoType = IBInfo;
    using GridType = IBGrid;
    using MetaType = CylinderData;
    using DataType = IBDataHolder<Cylinder>;

    static const std::string identifier() { return "Cylinder"; }
};

} // namespace ib
} // namespace amr_wind

#endif /* CYLINDER_H */
