#ifndef Joukowsky_H_
#define Joukowsky_H_

#include "amr-wind/wind_energy/actuator/disk/disk_types.H"
#include "amr-wind/wind_energy/actuator/disk/ActuatorDisk.H"
#include "amr-wind/core/vs/vector_space.H"

namespace amr_wind {
namespace actuator {

struct JoukowskyData : public DiskBaseData
{
    RealList angular_velocity;
    RealList tip_correction;
    RealList root_correction;
    amrex::Real current_angular_velocity{0.0};
};

struct Joukowsky : public DiskType
{
    using InfoType = ActInfo;
    using GridType = ActGrid;
    using MetaType = JoukowskyData;
    using DataType = ActDataHolder<Joukowsky>;

    static std::string identifier() { return "Joukowsky"; }
};

} // namespace actuator
} // namespace amr_wind

#endif /* Joukowsky_H_ */
