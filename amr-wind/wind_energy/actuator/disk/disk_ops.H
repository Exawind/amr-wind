#ifndef DISK_OPS_H_
#define DISK_OPS_H_

#include "amr-wind/wind_energy/actuator/disk/ActuatorDisk.H"
#include "amr-wind/wind_energy/actuator/actuator_ops.H"
#include "amr-wind/utilities/tensor_ops.H"

namespace amr_wind {
namespace actuator {
namespace disk {
// refactor these when we add more disk types later
void prepare_netcdf_file(
    const std::string&, const DiskBaseData&, const ActInfo&, const ActGrid&);

void write_netcdf(
    const std::string&,
    const DiskBaseData&,
    const ActInfo&,
    const ActGrid&,
    const amrex::Real);
} // namespace disk

namespace ops {

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

void check_for_parse_conflicts(const utils::ActParser& pp);

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
    DiskBaseData& meta,
    VecList& points,
    const vs::Vector& cylAxis,
    const int offset,
    const double dOffset);

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
        meta.disk_velocity = {0.0, 0.0, 0.0};
        auto& refVel = meta.reference_velocity;
        auto& diskVel = meta.disk_velocity;
        const int np = meta.num_vel_pts / 2;
        for (int i = 0; i < np; i++) {
            refVel = refVel + grid.vel[i];
            diskVel = diskVel + grid.vel[i + np];
        }
        refVel /= np;
        diskVel /= np;
    }
};
} // namespace ops
} // namespace actuator
} // namespace amr_wind

#endif /* DISK_OPS_H_ */
