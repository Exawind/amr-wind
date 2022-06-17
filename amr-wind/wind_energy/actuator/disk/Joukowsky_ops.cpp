#include "amr-wind/wind_energy/actuator/disk/Joukowsky_ops.H"
#include "amr-wind/utilities/ncutils/nc_interface.H"
#include "amr-wind/utilities/io_utils.H"

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
    pp.get("num_blades", meta.num_blades);
    meta.num_force_pts = meta.num_vel_pts_r * meta.num_vel_pts_t;
    meta.num_vel_pts = meta.num_force_pts * 2;
    meta.dr = 0.5 * meta.diameter / meta.num_vel_pts_r;
}

void parse_and_gather_params(const utils::ActParser& pp, JoukowskyData& data)
{
    check_for_parse_conflicts(pp);
    optional_parameters(data, pp);
    required_parameters(data, pp);
    ops::base::final_checks(data);
}

void prepare_netcdf_file(
    const std::string& name,
    const JoukowskyData& data,
    const ActInfo& info,
    const ActGrid& grid)
{
#ifdef AMR_WIND_USE_NETCDF
    using dvec = std::vector<double>;
    // Only root process handles I/O
    if (info.root_proc != amrex::ParallelDescriptor::MyProc()) return;
    auto ncf = ncutils::NCFile::create(name, NC_CLOBBER | NC_NETCDF4);
    const std::string nt_name = "num_time_steps";
    const std::string np_name = "num_actuator_points";
    const std::string nv_name = "num_velocity_points";

    ncf.enter_def_mode();
    ncf.put_attr("title", "AMR-Wind ActuatorDisk actuator output");
    ncf.put_attr("version", ioutils::amr_wind_version());
    ncf.put_attr("created_on", ioutils::timestamp());
    ncf.def_dim(nt_name, NC_UNLIMITED);
    ncf.def_dim("ndim", AMREX_SPACEDIM);

    auto grp = ncf.def_group(info.label);
    // clang-format off
    grp.put_attr("normal",
        dvec{data.normal_vec.x(), data.normal_vec.y(), data.normal_vec.z()});
    grp.put_attr("sample_normal",
        dvec{data.sample_vec.x(), data.sample_vec.y(), data.sample_vec.z()});
    // clang-format on
    grp.put_attr("diameter", dvec{data.diameter});
    grp.put_attr("epsilon", dvec{data.epsilon});
    grp.put_attr("sample_diameters", dvec{data.diameters_to_sample});
    grp.def_dim(np_name, data.num_force_pts);
    grp.def_dim(nv_name, data.num_vel_pts);
    grp.def_var("time", NC_DOUBLE, {nt_name});
    grp.def_var("xyz", NC_DOUBLE, {np_name, "ndim"});
    grp.def_var("xyz_v", NC_DOUBLE, {nv_name, "ndim"});
    grp.def_var("vref", NC_DOUBLE, {nt_name, "ndim"});
    grp.def_var("vdisk", NC_DOUBLE, {nt_name, "ndim"});
    grp.def_var("ct", NC_DOUBLE, {nt_name});
    grp.def_var("cp", NC_DOUBLE, {nt_name});
    grp.def_var("density", NC_DOUBLE, {nt_name});
    grp.def_var("total_disk_force", NC_DOUBLE, {nt_name, "ndim"});
    grp.def_var("angular_velocity", NC_DOUBLE, {nt_name});
    ncf.exit_def_mode();

    {
        {
            const size_t npts = static_cast<size_t>(data.num_force_pts);
            const std::vector<size_t> start{0, 0};
            const std::vector<size_t> count{npts, AMREX_SPACEDIM};
            grp.var("xyz").put(&(grid.pos[0][0]), start, count);
        }
        {
            const size_t npts = static_cast<size_t>(data.num_vel_pts);
            const std::vector<size_t> start{0, 0};
            const std::vector<size_t> count{npts, AMREX_SPACEDIM};
            grp.var("xyz_v").put(&(grid.vel_pos[0][0]), start, count);
        }
    }
#else
    amrex::ignore_unused(name, data, info, grid);
#endif
}

void write_netcdf(
    const std::string& name,
    const JoukowskyData& data,
    const ActInfo& info,
    const ActGrid& /*unused*/,
    const amrex::Real time)
{
#ifdef AMR_WIND_USE_NETCDF
    // Only root process handles I/O
    if (info.root_proc != amrex::ParallelDescriptor::MyProc()) return;
    auto ncf = ncutils::NCFile::open(name, NC_WRITE);
    const std::string nt_name = "num_time_steps";
    // Index of next timestep
    const size_t nt = ncf.dim(nt_name).len();
    auto grp = ncf.group(info.label);
    grp.var("time").put(&time, {nt}, {1});
    grp.var("vref").put(
        &(data.reference_velocity[0]), {nt, 0}, {1, AMREX_SPACEDIM});
    grp.var("vdisk").put(
        &(data.mean_disk_velocity[0]), {nt, 0}, {1, AMREX_SPACEDIM});
    grp.var("ct").put(&data.current_ct, {nt}, {1});
    grp.var("cp").put(&data.current_cp, {nt}, {1});
    grp.var("density").put(&data.density, {nt}, {1});
    grp.var("total_disk_force")
        .put(&data.disk_force[0], {nt}, {1, AMREX_SPACEDIM});
    grp.var("angular_velocity").put(&data.current_angular_velocity, {nt}, {1});
#else
    amrex::ignore_unused(name, data, info, time);
#endif
}

} // namespace joukowsky
} // namespace ops
} // namespace actuator
} // namespace amr_wind