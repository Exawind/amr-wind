#include "amr-wind/wind_energy/actuator/disk/Joukowsky_ops.H"

namespace amr_wind {
namespace actuator {
namespace ops {
namespace joukowsky {
void check_for_parse_conflicts(const utils::ActParser& pp)
{
    auto error_collector = ops::base::check_for_parse_conflicts(pp);
}
void optional_parameters(JoukowskyData& meta, const utils::ActParser& pp)
{
    ops::base::optional_parameters(meta, pp);
    pp.queryarr("angular_velocity", meta.angular_velocity);
    assert(!meta.angular_velocity.empty());
}
void required_parameters(JoukowskyData& meta, const utils::ActParser& pp)
{
    ops::base::required_parameters(meta, pp);
    pp.get("num_points_t", meta.num_vel_pts_t);
    pp.get("num_points_r", meta.num_vel_pts_r);
    meta.num_force_pts = meta.num_vel_pts_r * meta.num_vel_pts_t;
    meta.num_vel_pts = meta.num_force_pts * 2;
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