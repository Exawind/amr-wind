#ifndef CONSTANTCT_H_
#define CONSTANTCT_H_

#include "amr-wind/wind_energy/actuator/disk/disk_types.H"
#include "amr-wind/core/vs/vector_space.H"

namespace amr_wind {
namespace actuator {

struct ConstantCtData
{
    int num_points;
    vs::Vector center;
    vs::Vector normal;
    amrex::Real epsilon;
    amrex::Real diameter;
    amrex::Real dr;
};

struct ConstantCt : public DiskType
{
    using InfoType = ActInfo;
    using GridType = ActGrid;
    using MetaType = ConstantCtData;
    using DataType = ActDataHolder<ConstantCt>;

    static const std::string identifier() { return "ConstantCt"; }
};

} // namespace actuator
} // namespace amr_wind

#endif /* CONSTANTCT_H_ */
