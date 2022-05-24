#include "amr-wind/wind_energy/actuator/disk/Joukowsky_ops.H"

namespace amr_wind {
namespace actuator {
namespace ops {
namespace joukowsky {
void check_for_parse_conflicts(const utils::ActParser& pp)
{
    auto error_collector = ops::base::check_for_parse_conflicts(pp);
    // clang-format off
    base::collect_parse_dependencies(pp, "thrust_coeff", "angular_velocity", error_collector);
    base::collect_parse_dependencies(pp, "use_root_correction", "vortex_core_size", error_collector);
    // clang-format on
}
void optional_parameters(JoukowskyData& meta, const utils::ActParser& pp)
{
    ops::base::optional_parameters(meta, pp);
    pp.query("use_root_correction", meta.use_root_correction);
    pp.query("use_tip_correction", meta.use_tip_correction);
    pp.query("vortex_core_size", meta.vortex_core_size);
    pp.query("root_correction_exponent", meta.root_correction_exponent);
    pp.query("root_correction_coefficient", meta.root_correction_coefficient);
}
void required_parameters(JoukowskyData& meta, const utils::ActParser& pp)
{
    ops::base::required_parameters(meta, pp);
    pp.get("num_points_t", meta.num_vel_pts_t);
    pp.get("num_points_r", meta.num_vel_pts_r);
    pp.getarr("angular_velocity", meta.angular_velocity);
    meta.num_force_pts = meta.num_vel_pts_r * meta.num_vel_pts_t;
    meta.num_vel_pts = meta.num_force_pts * 2;
    meta.dr = 0.5 * meta.diameter / meta.num_vel_pts_r;
}
void parse_and_gather_params(const utils::ActParser& pp, JoukowskyData& meta)
{
    check_for_parse_conflicts(pp);
    optional_parameters(meta, pp);
    required_parameters(meta, pp);
    ops::base::final_checks(meta);
}
} // namespace joukowsky
} // namespace ops
} // namespace actuator
} // namespace amr_wind