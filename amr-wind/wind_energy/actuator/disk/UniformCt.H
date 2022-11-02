#ifndef UNIFORMCT_H_
#define UNIFORMCT_H_

#include "amr-wind/wind_energy/actuator/disk/ActuatorDisk.H"

namespace amr_wind::actuator {

struct UniformCtData : public DiskBaseData
{};

struct UniformCt : public DiskType
{
    using InfoType = ActInfo;
    using GridType = ActGrid;
    using MetaType = UniformCtData;
    using DataType = ActDataHolder<UniformCt>;

    static std::string identifier() { return "UniformCt"; }
};

} // namespace amr_wind::actuator

#endif /* UNIFORMCT_H_ */
