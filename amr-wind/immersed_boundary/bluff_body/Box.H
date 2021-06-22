#ifndef BOX_H
#define BOX_H

#include "amr-wind/immersed_boundary/bluff_body/BluffBody.H"

namespace amr_wind {
namespace ib {

struct BoxData : public BluffBodyBaseData
{
    vs::Vector center_loc{0.0, 0.0, 0.0};
    amrex::Real length{0.1};
    amrex::Real width{0.1};
    amrex::Real height{0.1};

    // force acting on immersed boundary
    amrex::Real fx{0.0};
    amrex::Real fy{0.0};
    amrex::Real fz{0.0};
};

struct Box : public BluffBodyType
{
    using InfoType = IBInfo;
    using MetaType = BoxData;
    using DataType = IBDataHolder<Box>;

    static const std::string identifier() { return "Box"; }
};

} // namespace ib
} // namespace amr_wind

#endif /* BOX_H */
