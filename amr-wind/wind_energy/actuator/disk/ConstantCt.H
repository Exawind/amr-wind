#ifndef CONSTANTCT_H_
#define CONSTANTCT_H_

#include "amr-wind/wind_energy/actuator/disk/disk_types.H"
#include "amr-wind/core/vs/vector_space.H"

namespace amr_wind {
namespace actuator {

struct ConstantCtData
{
    int num_force_pts;
    int num_vel_pts_r;
    int num_vel_pts_t;
    int num_vel_pts;
    vs::Vector center;
    vs::Vector normal_vec;
    vs::Vector coplanar_vec;
    vs::Vector reference_velocity;
    amrex::Real epsilon;
    amrex::Real diameter;
    amrex::Real dr;
    amrex::Real density;
    amrex::Real thrust_coeff;
    amrex::Real diameters_to_sample;
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
