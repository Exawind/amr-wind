#include "amr-wind/wind_energy/actuator/disk/Joukowski_ops.H"

namespace amr_wind {
namespace actuator {
namespace ops {
namespace joukowski {
void check_for_parse_conflicts(const utils::ActParser& pp)
{
    auto error_collector = ops::base::check_for_parse_conflicts(pp);
}
void optional_parameters(JoukowskiData& meta, const utils::ActParser& pp)
{
    ops::base::optional_parameters(meta, pp);
}
void required_parameters(JoukowskiData& meta, const utils::ActParser& pp)
{
    ops::base::required_parameters(meta, pp);
    pp.get("num_points_t", meta.num_vel_pts_t);
    pp.get("num_points_r", meta.num_vel_pts_r);
    meta.num_force_pts = meta.num_vel_pts_r * meta.num_vel_pts_t;
}
void parse_and_gather_params(const utils::ActParser& pp, JoukowskiData& meta)
{
    check_for_parse_conflicts(pp);
    optional_parameters(meta, pp);
    required_parameters(meta, pp);
    ops::base::final_checks(meta);
}
} // namespace joukowski
} // namespace ops
} // namespace actuator
} // namespace amr_wind