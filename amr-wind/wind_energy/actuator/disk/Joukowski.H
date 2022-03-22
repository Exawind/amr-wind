#ifndef JOUKOWSKI_H_
#define JOUKOWSKI_H_

#include "amr-wind/wind_energy/actuator/disk/disk_types.H"
#include "amr-wind/wind_energy/actuator/disk/ActuatorDisk.H"
#include "amr-wind/core/vs/vector_space.H"

namespace amr_wind {
namespace actuator {

struct JoukowskiData : public DiskBaseData
{
    amrex::Real power_coeff;
    amrex::Real power;
    amrex::Real moment;
    amrex::Real omega;
    amrex::Real tip_speed_ratio;
    RealList tip_correction;
    RealList root_correction;
};

struct Joukowski : public DiskType
{
    using InfoType = ActInfo;
    using GridType = ActGrid;
    using MetaType = JoukowskiData;
    using DataType = ActDataHolder<Joukowski>;

    static const std::string identifier() { return "Joukowski"; }
};

} // namespace actuator
} // namespace amr_wind

#endif /* JOUKOWSKI_H_ */
