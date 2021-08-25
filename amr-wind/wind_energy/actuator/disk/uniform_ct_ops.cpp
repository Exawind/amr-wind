#include "amr-wind/wind_energy/actuator/disk/uniform_ct_ops.H"
#include "amr-wind/utilities/ncutils/nc_interface.H"
#include "amr-wind/utilities/io_utils.H"
namespace amr_wind {
namespace actuator {
namespace disk {
void prepare_netcdf_file(
    const std::string& ncfile,
    const UniformCtData& meta,
    const ActInfo& info,
    const ActGrid& grid)
{
#ifdef AMR_WIND_USE_NETCDF
    using dvec = std::vector<double>;
    // Only root process handles I/O
    if (info.root_proc != amrex::ParallelDescriptor::MyProc()) return;
    auto ncf = ncutils::NCFile::create(ncfile, NC_CLOBBER | NC_NETCDF4);
    const std::string nt_name = "num_time_steps";
    const std::string np_name = "num_actuator_points";
    const std::string nv_name = "num_velocity_points";

    ncf.enter_def_mode();
    ncf.put_attr("title", "AMR-Wind UniformCtDisk actuator output");
    ncf.put_attr("version", ioutils::amr_wind_version());
    ncf.put_attr("created_on", ioutils::timestamp());
    ncf.def_dim(nt_name, NC_UNLIMITED);
    ncf.def_dim("ndim", AMREX_SPACEDIM);

    auto grp = ncf.def_group(info.label);
    // clang-format off
    grp.put_attr("normal",
        dvec{meta.normal_vec.x(), meta.normal_vec.y(), meta.normal_vec.z()});
    grp.put_attr("sample_normal",
        dvec{meta.sample_vec.x(), meta.sample_vec.y(), meta.sample_vec.z()});
    // clang-format on
    grp.put_attr("diameter", dvec{meta.diameter});
    grp.put_attr("epsilon", dvec{meta.epsilon});
    grp.put_attr("sample_diameters", dvec{meta.diameters_to_sample});
    grp.def_dim(np_name, meta.num_force_pts);
    grp.def_dim(nv_name, meta.num_vel_pts);
    grp.def_var("time", NC_DOUBLE, {nt_name});
    grp.def_var("xyz", NC_DOUBLE, {np_name, "ndim"});
    grp.def_var("xyz_v", NC_DOUBLE, {nv_name, "ndim"});
    grp.def_var("vref", NC_DOUBLE, {nt_name, "ndim"});
    grp.def_var("vdisk", NC_DOUBLE, {nt_name, "ndim"});
    grp.def_var("ct", NC_DOUBLE, {nt_name});
    grp.def_var("density", NC_DOUBLE, {nt_name});
    ncf.exit_def_mode();

    {
        {
            const size_t npts = static_cast<size_t>(meta.num_force_pts);
            const std::vector<size_t> start{0, 0};
            const std::vector<size_t> count{npts, AMREX_SPACEDIM};
            grp.var("xyz").put(&(grid.pos[0][0]), start, count);
        }
        {
            const size_t npts = static_cast<size_t>(meta.num_vel_pts);
            const std::vector<size_t> start{0, 0};
            const std::vector<size_t> count{npts, AMREX_SPACEDIM};
            grp.var("xyz_v").put(&(grid.vel_pos[0][0]), start, count);
        }
    }
#else
    amrex::ignore_unused(ncfile, meta, info, grid);
#endif
}

void write_netcdf(
    const std::string& ncfile,
    const UniformCtData& meta,
    const ActInfo& info,
    const ActGrid&,
    const amrex::Real time)
{
#ifdef AMR_WIND_USE_NETCDF
    // Only root process handles I/O
    if (info.root_proc != amrex::ParallelDescriptor::MyProc()) return;
    auto ncf = ncutils::NCFile::open(ncfile, NC_WRITE);
    const std::string nt_name = "num_time_steps";
    // Index of next timestep
    const size_t nt = ncf.dim(nt_name).len();
    auto grp = ncf.group(info.label);
    grp.var("time").put(&time, {nt}, {1});
    grp.var("vref").put(
        &(meta.reference_velocity[0]), {nt, 0}, {1, AMREX_SPACEDIM});
    grp.var("vdisk").put(
        &(meta.disk_velocity[0]), {nt, 0}, {1, AMREX_SPACEDIM});
    grp.var("ct").put(&meta.current_ct, {nt}, {1});
    grp.var("density").put(&meta.density, {nt}, {1});
#else
    amrex::ignore_unused(ncfile, meta, info, time);
#endif
}
} // namespace disk
namespace ops {

vs::Vector get_east_orientation()
{
    utils::ActParser pp("Coriolis.Forcing", "Coriolis");
    vs::Vector east;
    if (pp.contains("east_vector")) {
        pp.get("east_vector", east);
    } else {
        east = vs::Vector::ihat();
    }
    return east;
}

vs::Vector get_north_orientation()
{
    utils::ActParser pp("Coriolis.Forcing", "Coriolis");
    vs::Vector north;
    if (pp.contains("north_vector")) {
        pp.get("north_vector", north);
    } else {
        north = vs::Vector::jhat();
    }
    return north;
}

void collect_parse_conflicts(
    const utils::ActParser& pp,
    const std::string& p1,
    const std::string& p2,
    std::ostringstream& ss)
{
    if (pp.contains(p1) && pp.contains(p2)) {
        ss << "UniformCt Conflict: " << p1 << " and " << p2 << std::endl;
    }
}
void collect_parse_dependencies(
    const utils::ActParser& pp,
    const std::string& p1,
    const std::string& p2,
    std::ostringstream& ss)
{
    if (pp.contains(p1) && !pp.contains(p2))
        ss << "UniformCt Dependency Missing: " << p1 << " and " << p2
           << std::endl;
    if (!pp.contains(p1) && pp.contains(p2))
        ss << "UniformCt Dependency Missing: " << p1 << " and " << p2
           << std::endl;
}

void required_parameters(UniformCt::MetaType& meta, const utils::ActParser& pp)
{
    pp.get("num_force_points", meta.num_force_pts);
    pp.get("epsilon", meta.epsilon);
    pp.get("rotor_diameter", meta.diameter);
    pp.getarr("thrust_coeff", meta.thrust_coeff);
}

void optional_parameters(UniformCt::MetaType& meta, const utils::ActParser& pp)
{
    // no logic is required to check for conflicts in this function. all
    // conflicts should be specified in the check_for_parse_conflicts function
    // where we can check for all conflicts and dump them all out together.
    // hopefully this way users will only have to update mistakes in one
    // iteration

    meta.normal_vec = get_north_orientation();
    meta.sample_vec = get_north_orientation();

    if (pp.contains("base_position")) {
        amrex::Real hub;
        vs::Vector base;
        pp.get("base_position", base);
        pp.get("hub_height", hub);
        meta.center = base + (hub * vs::Vector::khat());
    }
    pp.query("disk_center", meta.center);
    pp.query("disk_normal", meta.normal_vec);
    pp.query("density", meta.density);
    pp.query("diameters_to_sample", meta.diameters_to_sample);

    // make sure we compute normal vec contribution from tilt before yaw
    // since we won't know a reference axis to rotate for tilt after
    // yawing
    const auto east = get_east_orientation();

    auto normalRotOp = vs::Tensor::I();
    auto sampleRotOp = vs::Tensor::I();

    if (pp.contains("tilt")) {
        amrex::Real tilt;
        pp.get("tilt", tilt);
        const auto tilt_op = vs::quaternion(east, -tilt);
        normalRotOp = normalRotOp & tilt_op;
    }

    if (pp.contains("yaw")) {
        amrex::Real yaw;
        pp.get("yaw", yaw);
        // use negative yaw angle to match convention of clockwise rotation
        normalRotOp = normalRotOp & vs::zrot(-yaw);
    }

    bool user_specified_sample_vec = false;

    pp.query("sample_normal", meta.sample_vec);

    if (pp.contains("sample_tilt")) {
        user_specified_sample_vec = true;
        amrex::Real tilt;
        pp.get("sample_tilt", tilt);
        const auto tilt_op = vs::quaternion(east, -tilt);
        sampleRotOp = sampleRotOp & tilt_op;
    }
    if (pp.contains("sample_yaw")) {
        user_specified_sample_vec = true;
        amrex::Real yaw;
        pp.get("sample_yaw", yaw);
        // use negative yaw angle to match convention of clockwise rotation
        sampleRotOp = sampleRotOp & vs::zrot(-yaw);
    }

    meta.normal_vec = meta.normal_vec & normalRotOp;
    meta.sample_vec = meta.sample_vec & sampleRotOp;
    if (!pp.contains("sample_normal") && !user_specified_sample_vec) {
        // if none of the sample vec operations are supplied by user
        // then just use normal this will be overwritten if the user wants to
        // specify it themselves
        meta.sample_vec = meta.normal_vec;
    }

    // params for the sampling disk
    pp.query("num_vel_points_r", meta.num_vel_pts_r);
    pp.query("num_vel_points_t", meta.num_vel_pts_t);
    // 2x 1 for sampling up stream and one for sampling at the disk
    meta.num_vel_pts = meta.num_vel_pts_r * meta.num_vel_pts_t * 2;

    // ensure any computed vectors are normalized
    meta.normal_vec.normalize();
    meta.sample_vec.normalize();

    pp.queryarr("wind_speed", meta.table_velocity);
    if (meta.thrust_coeff.size() == 1) {
        meta.current_ct = meta.thrust_coeff[0];
        meta.table_velocity = {1.0};
    }
}

void check_for_parse_conflicts(const utils::ActParser& pp)
{
    std::ostringstream error_collector;

    // clang-format off
    collect_parse_conflicts(pp, "disk_normal", "yaw", error_collector);
    collect_parse_conflicts(pp, "disk_normal", "tilt", error_collector);
    collect_parse_conflicts(pp, "sample_normal", "sample_yaw", error_collector);
    collect_parse_conflicts(pp, "sample_normal", "sample_tilt", error_collector);
    collect_parse_conflicts(pp, "disk_center", "base_position", error_collector);
    collect_parse_conflicts(pp, "disk_center", "hub_height", error_collector);
    collect_parse_dependencies(pp, "base_position", "hub_height", error_collector);
    // clang-format on
    RealList ct;
    pp.getarr("thrust_coeff", ct);
    if (ct.size() > 1) {
        if (!pp.contains("wind_speed")) {
            error_collector
                << "UniformCt Dependency Missing: wind_speed is required when "
                   "there is more than 1 entry for thrust_coeff"
                << std::endl;
        }
    }

    if (pp.contains("wind_speed")) {
        RealList vel;
        pp.getarr("wind_speed", vel);
        if (vel.size() != ct.size())
            error_collector << "UniformCt Conflict: wind_speed and "
                               "thrust_coeff must have the same number of "
                               "values"
                            << std::endl;
    }

    if (!error_collector.str().empty())
        amrex::Abort(
            "Errors found while parsing in Actuator.UniformCt:\n" +
            error_collector.str());
}

void compute_and_normalize_coplanar_vector(UniformCt::MetaType& meta)
{
    const amrex::Real radius = meta.diameter * 0.5;
    meta.dr = radius / meta.num_force_pts;

    // ensure normal is normalized
    meta.normal_vec.normalize();
    const vs::Vector& norm = meta.normal_vec;

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        vs::mag_sqr(vs::Vector::khat() ^ norm) >
            vs::DTraits<amrex::Real>::eps(),
        "A constant Ct disk normal is too close to vertical.");

    // compute a coplanar vector that resides in the same plane as the disk
    // we will use this vector for the bounding box and the force point
    // locations
    meta.coplanar_vec = norm ^ vs::Vector::khat();
    meta.coplanar_vec.normalize();
    meta.sample_vec.normalize();
}

void final_checks(const UniformCt::MetaType& meta)
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        meta.num_vel_pts > 0,
        "num_vel_points_r and num_vel_points_t must both be >=1");
}

amrex::RealBox compute_bounding_box(const UniformCt::MetaType& meta)
{
    auto& norm = meta.normal_vec;
    auto& cVec = meta.coplanar_vec;

    const auto& cc = meta.center;
    const amrex::Real nl = meta.epsilon * 3.0; // length scale in normal dir
    const amrex::Real dl =
        meta.diameter * 0.5 + meta.dr * 2.0; // length scale in plane of disk
    const auto dvec = norm * nl + cVec * dl + vs::Vector::khat() * dl;
    const auto p1 = cc - dvec; // front
    const auto p2 = cc + dvec; // back
    return amrex::RealBox(
        amrex::min(p1.x(), p2.x()), amrex::min(p1.y(), p2.y()),
        amrex::min(p1.z(), p2.z()), amrex::max(p1.x(), p2.x()),
        amrex::max(p1.y(), p2.y()), amrex::max(p1.z(), p2.z()));
}

void do_parse_based_computations(UniformCt::DataType& data)
{
    auto& meta = data.meta();
    auto& info = data.info();

    compute_and_normalize_coplanar_vector(meta);
    info.bound_box = compute_bounding_box(meta);
}

void compute_disk_points(
    UniformCt::MetaType& meta,
    VecList& points,
    const vs::Vector& cylAxis,
    const int offset,
    const double dOffset)
{
    const auto& cc = meta.center;
    // we define points as if the disk faces the standard cyindrical coordinates
    // and then rotate it into position

    // first get an angle orthogonal to centerline of the cylinder and the z
    // direction (standard centerline for a cylinder). this will be the vector
    // we will rotate about
    auto rotVec = vs::Vector::khat() ^ cylAxis;
    rotVec.normalize();

    // next compute the angle between the center vec and z axis
    const amrex::Real angle =
        ::amr_wind::utils::degrees(std::acos((cylAxis & vs::Vector::khat())));

    // finally get the rotation matrix given that angle and rotation vector
    const auto rotMatrix = vs::quaternion(rotVec, angle);

    const vs::VectorT<int> nvp = {meta.num_vel_pts_r, meta.num_vel_pts_t, 1};

    const amrex::Real dr = meta.diameter / nvp.x();
    const amrex::Real dt = ::amr_wind::utils::two_pi() / nvp.y();
    const amrex::Real du = dOffset * meta.diameter;

    int ip = offset;
    for (int i = 0; i < nvp.x(); i++) {
        const amrex::Real r = dr * i;
        for (int j = 0; j < nvp.y(); j++, ip++) {
            const amrex::Real theta = j * dt;
            vs::Vector refPoint = {
                r * std::cos(theta), r * std::sin(theta), du};
            points[ip] = (refPoint & rotMatrix) + cc;
        }
    }
}
} // namespace ops
} // namespace actuator
} // namespace amr_wind
