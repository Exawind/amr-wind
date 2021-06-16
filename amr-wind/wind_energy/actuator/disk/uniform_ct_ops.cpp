#include "amr-wind/wind_energy/actuator/disk/uniform_ct_ops.H"
namespace amr_wind {
namespace actuator {
namespace ops {
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
    if (!pp.contains(p1) && !pp.contains(p2))
        ss << "UniformCt Dependency Missing: " << p1 << " and " << p2
           << std::endl;
}

void required_parameters(UniformCt::MetaType& meta, const utils::ActParser& pp)
{
    pp.get("num_force_points", meta.num_force_pts);
    {
        amrex::Real hub;
        vs::Vector base;
        pp.get("base_position", base);
        pp.get("hub_height", hub);
        meta.center = base + (hub * vs::Vector::khat());
    }
    pp.get("epsilon", meta.epsilon);
    pp.get("rotor_diameter", meta.diameter);
    pp.get("thrust_coeff", meta.thrust_coeff);
}

void optional_parameters(UniformCt::MetaType& meta, const utils::ActParser& pp)
{
    pp.query("disk_normal", meta.normal_vec);
    pp.query("density", meta.density);
    pp.query("diameters_to_sample", meta.diameters_to_sample);

    if (pp.contains("yaw")) {
        amrex::Real yaw;
        pp.get("yaw", yaw);
        meta.normal_vec = vs::Vector::ihat() & vs::zrot(yaw);
    }
    if (pp.contains("sample_normal")) {
        pp.get("sample_normal", meta.sample_vec);
    } else if (pp.contains("sample_yaw")) {
        amrex::Real sYaw;
        pp.get("sample_yaw", sYaw);
        meta.sample_vec = vs::Vector::ihat() & vs::zrot(sYaw);
    } else {
        meta.sample_vec = meta.normal_vec;
    }
    pp.query("num_vel_points_r", meta.num_vel_pts_r);
    pp.query("num_vel_points_t", meta.num_vel_pts_t);
    meta.num_vel_pts = meta.num_vel_pts_r * meta.num_vel_pts_t;

    // ensure any computed vectors are normalized
    meta.normal_vec.normalize();
    meta.sample_vec.normalize();
}

void check_for_parse_conflicts(const utils::ActParser& pp)
{
    std::ostringstream error_collector;

    collect_parse_conflicts(pp, "disk_normal", "yaw", error_collector);
    collect_parse_conflicts(pp, "sample_normal", "sample_yaw", error_collector);

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
} // namespace ops
} // namespace actuator
} // namespace amr_wind