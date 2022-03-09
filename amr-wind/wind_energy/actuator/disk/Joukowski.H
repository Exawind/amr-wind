#ifndef JOUKOWSKI_H_
#define JOUKOWSKI_H_

#include "amr-wind/wind_energy/actuator/disk/disk_types.H"
#include "amr-wind/wind_energy/actuator/disk/ActuatorDisk.H"
#include "amr-wind/core/vs/vector_space.H"

namespace amr_wind {
namespace actuator {

struct JouwkoskiData : public DiskBaseData
{};

struct Joukowski : public DiskType
{
    using InfoType = ActInfo;
    using GridType = ActGrid;
    using MetaType = JouwkoskiData;
    using DataType = ActDataHolder<Joukowski>;

    static const std::string identifier() { return "Joukowski"; }
};

} // namespace actuator
} // namespace amr_wind

#endif /* JOUKOWSKI_H_ */
