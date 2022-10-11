#ifndef Joukowsky_H_
#define Joukowsky_H_

#include "amr-wind/wind_energy/actuator/disk/disk_types.H"
#include "amr-wind/wind_energy/actuator/disk/ActuatorDisk.H"
#include "amr-wind/core/vs/vector_space.H"

namespace amr_wind::actuator {

struct JoukowskyData : public DiskBaseData
{
    RealList angular_velocity;
    RealList tip_correction;
    RealList root_correction;
    RealList f_normal;
    RealList f_theta;
    amrex::Real current_angular_velocity{0.0};
    amrex::Real current_tip_speed_ratio{0.0};
    amrex::Real vortex_core_size;
    amrex::Real current_cp;
    amrex::Real current_power;
    // --- Sorenson 2020 equation 10 constants ----
    amrex::Real root_correction_coefficient{2.335};
    amrex::Real root_correction_exponent{4.0};
    // --------------------------------------------
    int num_blades{3};
    vs::Vector disk_force{0.0, 0.0, 0.0};
    bool use_tip_correction{true};
    bool use_root_correction{true};
};

struct Joukowsky : public DiskType
{
    using InfoType = ActInfo;
    using GridType = ActGrid;
    using MetaType = JoukowskyData;
    using DataType = ActDataHolder<Joukowsky>;

    static std::string identifier() { return "Joukowsky"; }
};

} // namespace amr_wind::actuator

#endif /* Joukowsky_H_ */
