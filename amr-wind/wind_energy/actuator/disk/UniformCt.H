#ifndef UNIFORMCT_H_
#define UNIFORMCT_H_

#include "amr-wind/wind_energy/actuator/disk/disk_types.H"
#include "amr-wind/core/vs/vector_space.H"

namespace amr_wind {
namespace actuator {

struct UniformCtData
{
    int num_force_pts;
    int num_vel_pts{1};
    int num_vel_pts_r{1};
    int num_vel_pts_t{1};
    vs::Vector center;
    vs::Vector normal_vec{0, 1, 0};
    vs::Vector coplanar_vec{1, 0, 0};
    vs::Vector sample_vec{0, 1, 0};
    vs::Vector reference_velocity;
    vs::Vector disk_velocity;
    amrex::Real epsilon;
    amrex::Real diameter;
    amrex::Real dr;
    amrex::Real density;
    amrex::Real current_ct;
    amrex::Real diameters_to_sample{0.0};
    RealList thrust_coeff;
    RealList table_velocity;
};

struct UniformCt : public DiskType
{
    using InfoType = ActInfo;
    using GridType = ActGrid;
    using MetaType = UniformCtData;
    using DataType = ActDataHolder<UniformCt>;

    static const std::string identifier() { return "UniformCt"; }
};

} // namespace actuator
} // namespace amr_wind

#endif /* UNIFORMCT_H_ */
