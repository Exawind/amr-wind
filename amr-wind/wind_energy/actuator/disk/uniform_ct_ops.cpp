#include "amr-wind/wind_energy/actuator/disk/uniform_ct_ops.H"

namespace amr_wind {
namespace actuator {
namespace ops {
namespace uniformct {
void check_for_parse_conflicts(const utils::ActParser& pp)
{
    auto error_collector = ops::base::check_for_parse_conflicts(pp);
    ops::base::collect_parse_dependencies(
        pp, "num_theta_force_points", "spreading", error_collector);
}
void optional_parameters(DiskBaseData& meta, const utils::ActParser& pp)
{
    ops::base::optional_parameters(meta, pp);
    // params for the sampling disk
    pp.query("num_vel_points_r", meta.num_vel_pts_r);
    pp.query("num_vel_points_t", meta.num_vel_pts_t);
    // 2x 1 for sampling up stream and one for sampling at the disk
    meta.num_vel_pts = meta.num_vel_pts_r * meta.num_vel_pts_t * 2;
    pp.query("spreading_type", meta.spreading_type);
}
void required_parameters(DiskBaseData& meta, const utils::ActParser& pp)
{
    ops::base::required_parameters(meta, pp);
    pp.get("num_force_points", meta.num_force_pts);
    meta.dr = 0.5 * meta.diameter / meta.num_force_pts;
}
void parse_and_gather_params(const utils::ActParser& pp, UniformCtData& data)
{
    check_for_parse_conflicts(pp);
    optional_parameters(data, pp);
    required_parameters(data, pp);
    ops::base::final_checks(data);
}
} // namespace uniformct
} // namespace ops
} // namespace actuator
} // namespace amr_wind