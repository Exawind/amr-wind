#ifndef RECTANGULAR_PRISM_H
#define RECTANGULAR_PRISM_H

#include "amr-wind/immersed_boundary/geometry/BluffBody.H"

namespace amr_wind {
namespace ib {

struct RectangularPrismData : public BluffBodyBaseData
{
    vs::Vector center_loc{0.0, 0.0, 0.0};
    amrex::Real length{0.1};
    amrex::Real width{0.1};
    amrex::Real height{0.1};
};

struct RectangularPrism : public BluffBodyType
{
    using InfoType = IBInfo;
    using GridType = IBGrid;
    using MetaType = RectangularPrismData;
    using DataType = IBDataHolder<RectangularPrism>;

    static const std::string identifier() { return "RectangularPrism"; }
};

} // namespace ib
} // namespace amr_wind

#endif /* RECTANGULAR_PRISM_H */
