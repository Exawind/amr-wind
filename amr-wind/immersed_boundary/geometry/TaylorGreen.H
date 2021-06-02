#ifndef TAYLORGREEN_H
#define TAYLORGREEN_H

#include "amr-wind/immersed_boundary/geometry/BluffBody.H"

namespace amr_wind {
namespace ib {

struct TaylorGreenData : public BluffBodyBaseData
{
    vs::Vector center_loc{0.0, 0.0, 0.0};
    amrex::Real radius{0.1};
};

struct TaylorGreen : public BluffBodyType
{
    using InfoType = IBInfo;
    using GridType = IBGrid;
    using MetaType = TaylorGreenData;
    using DataType = IBDataHolder<TaylorGreen>;

    static const std::string identifier() { return "TaylorGreen"; }
};

} // namespace ib
} // namespace amr_wind

#endif /* TAYLORGREEN_H */
