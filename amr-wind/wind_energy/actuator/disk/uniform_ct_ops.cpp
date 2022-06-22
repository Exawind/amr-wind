#include "amr-wind/wind_energy/actuator/disk/uniform_ct_ops.H"
#include "amr-wind/wind_energy/actuator/ActParser.H"
#include <ostream>

namespace amr_wind {
namespace actuator {
namespace ops {
namespace uniformct {

void check_for_removed_syntax(
    const utils::ActParser& pp, std::ostringstream& stream)
{
    const std::string error_start = "ERROR ActuatorDisk:: ";
    if (pp.contains("num_vel_points_r")) {
        stream << error_start
               << "'num_vel_points_r' has been replaced with 'num_points_r'\n";
    }
    if (pp.contains("num_vel_points_t")) {
        stream << error_start
               << "'num_vel_points_t' has been replaced with 'num_points_t'\n";
    }
    if (pp.contains("num_force_points")) {
        stream << error_start
               << "'num_force_points' has been repalced with 'num_points_r'\n";
    }
    if (pp.contains("num_theta_force_points")) {
        stream << error_start
               << "'num_theta_force_points' has been replaced with "
                  "'num_points_t'\n";
    }
}

void check_for_parse_conflicts(const utils::ActParser& pp)
{
    auto error_collector = ops::base::check_for_parse_conflicts(pp);
    check_for_removed_syntax(pp, error_collector);
    ops::base::collect_parse_dependencies_one_way(
        pp, "num_points_t", "spreading", error_collector);
    ops::base::check_error_stream(error_collector);
}
void optional_parameters(DiskBaseData& meta, const utils::ActParser& pp)
{
    ops::base::optional_parameters(meta, pp);
    // params for the sampling disk
    pp.query("num_points_t", meta.num_vel_pts_t);
    // 2x 1 for sampling up stream and one for sampling at the disk
    meta.num_vel_pts = meta.num_vel_pts_r * meta.num_vel_pts_t * 2;
    pp.query("spreading_type", meta.spreading_type);
}
void required_parameters(DiskBaseData& meta, const utils::ActParser& pp)
{
    ops::base::required_parameters(meta, pp);
    pp.get("num_points_r", meta.num_force_pts);
    meta.num_vel_pts_r = meta.num_force_pts;
    meta.dr = 0.5 * meta.diameter / meta.num_force_pts;
}
void parse_and_gather_params(const utils::ActParser& pp, UniformCtData& data)
{
    check_for_parse_conflicts(pp);
    required_parameters(data, pp);
    optional_parameters(data, pp);
    ops::base::final_checks(data);
}
} // namespace uniformct
} // namespace ops
} // namespace actuator
} // namespace amr_wind