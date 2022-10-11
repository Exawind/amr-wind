#ifndef BOX_H
#define BOX_H

#include "amr-wind/immersed_boundary/bluff_body/BluffBody.H"

namespace amr_wind::ib {

struct BoxData : public BluffBodyBaseData
{
    vs::Vector center_loc{0.0, 0.0, 0.0};
    amrex::Real length{0.1};
    amrex::Real width{0.1};
    amrex::Real height{0.1};
};

struct Box : public BluffBodyType
{
    using InfoType = IBInfo;
    using MetaType = BoxData;
    using DataType = IBDataHolder<Box>;

    static std::string identifier() { return "Box"; }
};

} // namespace amr_wind::ib

#endif /* BOX_H */
