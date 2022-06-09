#ifndef DISK_OPS_H_
#define DISK_OPS_H_

#include "amr-wind/wind_energy/actuator/disk/ActuatorDisk.H"
#include "amr-wind/wind_energy/actuator/actuator_ops.H"
#include "amr-wind/utilities/tensor_ops.H"
#include "amr-wind/utilities/linear_interpolation.H"

namespace amr_wind {
namespace actuator {
namespace disk {
// refactor these when we add more disk types later
void prepare_netcdf_file(
    const std::string& name,
    const DiskBaseData& data,
    const ActInfo& info,
    const ActGrid& grid);

void write_netcdf(
    const std::string& name,
    const DiskBaseData& data,
    const ActInfo& info,
    const ActGrid& /*unused*/,
    const amrex::Real time);
} // namespace disk

namespace ops {
namespace base {
vs::Vector compute_coplanar_vector(const vs::Vector& normal);

void collect_parse_conflicts(
    const utils::ActParser& pp,
    const std::string& p1,
    const std::string& p2,
    std::ostringstream& ss);
void collect_parse_dependencies(
    const utils::ActParser& pp,
    const std::string& p1,
    const std::string& p2,
    std::ostringstream& ss);

void required_parameters(DiskBaseData& meta, const utils::ActParser& pp);

void optional_parameters(DiskBaseData& meta, const utils::ActParser& pp);

std::ostringstream check_for_parse_conflicts(const utils::ActParser& pp);

void compute_and_normalize_coplanar_vector(DiskBaseData& meta);

void final_checks(const DiskBaseData& meta);

amrex::RealBox compute_bounding_box(const DiskBaseData& meta);

template <typename T>
void do_parse_based_computations(ActDataHolder<T>& data)
{
    auto& meta = data.meta();
    auto& info = data.info();

    compute_and_normalize_coplanar_vector(meta);
    info.bound_box = compute_bounding_box(meta);
}

void compute_disk_points(
    const DiskBaseData& meta,
    VecList& points,
    const vs::Vector& cylAxis,
    const int offset,
    const double dOffset);

template <typename T>
void parse_and_gather_params(const utils::ActParser& pp, T& metaData)
{
    check_for_parse_conflicts(pp);
    required_parameters(metaData, pp);
    optional_parameters(metaData, pp);
}

inline amrex::Real compute_reference_velocity_sqr(DiskBaseData& data)
{
    const auto& normal = data.normal_vec;
    // only use velocity components normal to disk
    // make sure uInfSqr is always positive by squaring the projection of
    // reference_velocity onto the normal
    return (data.reference_velocity & normal) *
           (data.reference_velocity & normal);
}

// we use a setter function since we want to cache this variable for netcdf
// output
inline void
set_thrust_coefficient(DiskBaseData& data, const amrex::Real& uInfSqr)
{
    const amrex::Real uInfMag = std::sqrt(uInfSqr);
    data.current_ct = ::amr_wind::interp::linear(
        data.table_velocity, data.thrust_coeff, uInfMag);
}

} // namespace base

template <typename ActTrait>
struct UpdateVelOp<
    ActTrait,
    ActSrcDisk,
    std::enable_if_t<std::is_base_of<DiskType, ActTrait>::value>>
{
    void operator()(typename ActTrait::DataType& data)
    {
        auto& meta = data.meta();
        const auto& grid = data.grid();
        meta.reference_velocity = {0.0, 0.0, 0.0};
        meta.mean_disk_velocity = {0.0, 0.0, 0.0};
        auto& refVel = meta.reference_velocity;
        auto& diskVel = meta.mean_disk_velocity;
        const int np = meta.num_vel_pts / 2;
        for (int i = 0; i < np; i++) {
            refVel = refVel + grid.vel[i];
            diskVel = diskVel + grid.vel[i + np];
        }
        refVel /= np;
        diskVel /= np;
    }
};

template <typename ActTrait>
struct UpdatePosOp<
    ActTrait,
    ActSrcDisk,
    std::enable_if_t<std::is_base_of<DiskType, ActTrait>::value>>
{
    void operator()(typename ActTrait::DataType& /*data*/) {}
};
} // namespace ops
} // namespace actuator
} // namespace amr_wind

#endif /* DISK_OPS_H_ */
