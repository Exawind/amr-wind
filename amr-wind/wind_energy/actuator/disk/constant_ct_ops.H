#ifndef CONSTANT_CT_OPS_H_
#define CONSTANT_CT_OPS_H_

#include "amr-wind/wind_energy/actuator/disk/ConstantCt.H"
#include "amr-wind/wind_energy/actuator/actuator_ops.H"
#include "amr-wind/utilities/tensor_ops.H"

namespace amr_wind {
namespace actuator {
namespace ops {

namespace {
void collect_parse_conflicts(
    const utils::ActParser& pp,
    const std::string& p1,
    const std::string& p2,
    std::ostringstream& ss)
{
    if (pp.contains(p1) && pp.contains(p2)) {
        ss << "ConstantCt Conflict: " << p1 << " and " << p2 << std::endl;
    }
}
void collect_parse_dependencies(
    const utils::ActParser& pp,
    const std::string& p1,
    const std::string& p2,
    std::ostringstream& ss)
{
    if (pp.contains(p1) && !pp.contains(p2))
        ss << "ConstantCt Dependency Missing: " << p1 << " and " << p2
           << std::endl;
    if (!pp.contains(p1) && !pp.contains(p2))
        ss << "ConstantCt Dependency Missing: " << p1 << " and " << p2
           << std::endl;
}

void required_parameters(ConstantCt::MetaType& meta, const utils::ActParser& pp)
{
    pp.get("num_force_points", meta.num_force_pts);
    pp.get("disk_center", meta.center);
    pp.get("epsilon", meta.epsilon);
    pp.get("rotor_diameter", meta.diameter);
    pp.get("thrust_coeff", meta.thrust_coeff);
}

void optional_parameters(ConstantCt::MetaType& meta, const utils::ActParser& pp)
{
    pp.query("disk_normal", meta.normal_vec);
    pp.query("density", meta.density);
    pp.query("diameters_to_sample", meta.diameters_to_sample);

    if (pp.contains("yaw")) {
        amrex::Real yaw;
        pp.get("yaw", yaw);
        yaw = ::amr_wind::utils::radians(yaw);
        meta.normal_vec = vs::Vector::ihat() & vs::zrot(yaw);
    }
    if (pp.contains("sample_normal")) {
        pp.get("sample_normal", meta.sample_vec);
    } else {
        meta.sample_vec = meta.normal_vec;
    }
    pp.query("num_vel_points_r", meta.num_vel_pts_r);
    pp.query("num_vel_points_t", meta.num_vel_pts_t);
    meta.num_vel_pts = meta.num_vel_pts_r * meta.num_vel_pts_t;
}

void check_for_parse_conflicts(const utils::ActParser& pp)
{
    std::ostringstream error_collector;

    collect_parse_conflicts(pp, "disk_normal", "yaw", error_collector);

    if (!error_collector.str().empty())
        amrex::Abort(
            "Errors found while parsing a ConstantCt Actuator:\n" +
            error_collector.str());
}

void compute_and_normalize_vectors(ConstantCt::MetaType& meta)
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

void final_checks(const ConstantCt::MetaType& meta)
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        meta.num_vel_pts > 0,
        "num_vel_points_r and num_vel_points_t must both be >=1");
}

amrex::RealBox compute_bounding_box(const ConstantCt::MetaType& meta)
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

void do_parse_based_computations(ConstantCt::DataType& data)
{
    auto& meta = data.meta();
    auto& info = data.info();

    compute_and_normalize_vectors(meta);
    info.bound_box = compute_bounding_box(meta);
}
} // namespace
template <>
struct ReadInputsOp<ConstantCt, ActSrcDisk>
{
    void operator()(ConstantCt::DataType& data, const utils::ActParser& pp)
    {
        // TODO get density at probe location
        check_for_parse_conflicts(pp);
        auto& meta = data.meta();
        required_parameters(meta, pp);
        optional_parameters(meta, pp);
        final_checks(meta);
        do_parse_based_computations(data);
    }
};
template <>
struct InitDataOp<ConstantCt, ActSrcDisk>
{
    void operator()(ConstantCt::DataType& data)
    {
        auto& grid = data.grid();
        auto& meta = data.meta();

        // only resize the members we are going to use
        grid.pos.resize(meta.num_force_pts);
        grid.force.resize(meta.num_force_pts);
        grid.vel.resize(meta.num_vel_pts);
        grid.vel_pos.resize(meta.num_vel_pts);

        const auto& cVec = meta.coplanar_vec;
        const auto& sVec = meta.sample_vec;
        const auto& nVec = meta.normal_vec;

        const auto& cc = meta.center;
        {
            const auto& dr = meta.dr;
            for (int i = 0; i < meta.num_force_pts; ++i) {
                grid.pos[i] = cc + (i + 0.5) * dr * cVec;
            }
        }
        {
            const auto rotVec = vs::Vector::khat() ^ sVec;
            // vectors are normalized so magnitude of each is 1
            const amrex::Real angle = std::acos((sVec & vs::Vector::ihat()));
            const auto rotMatrix = vs::quaternion(rotVec, angle);

            const vs::Vector nvp = {meta.num_vel_pts_r, meta.num_vel_pts_t, 1};

            const amrex::Real dr = meta.diameter / nvp.x();
            const amrex::Real dt = ::amr_wind::utils::two_pi() / nvp.y();
            const amrex::Real du = meta.diameters_to_sample * meta.diameter;

            int ip = 0;
            for (int i = 0; i < nvp.x(); i++) {
                const amrex::Real r = dr * (i + 0.5);
                for (int j = 0; j < nvp.y(); j++, ip++) {
                    const amrex::Real theta = j * dt;
                    vs::Vector refPoint = {
                        r * std::cos(theta), r * std::sin(theta), du};
                    grid.vel_pos[ip] = (refPoint & rotMatrix) + cc;
                }
            }
        }
    }
};

template <>
struct UpdateVelOp<ConstantCt, ActSrcDisk>
{
    void operator()(ConstantCt::DataType& data)
    {
        auto& meta = data.meta();
        const auto& grid = data.grid();
        auto& refVel = meta.reference_velocity;
        for (auto&& v : grid.vel) {
            refVel = refVel + v;
        }
        refVel /= grid.vel.size();
    }
};

// compute the total force over a given radial section
// this will then be spread uniformly over that area in the
// projection step
template <>
struct ComputeForceOp<ConstantCt, ActSrcDisk>
{
    void operator()(ConstantCt::DataType& data)
    {
        auto& grid = data.grid();
        auto& ddata = data.meta();
        const amrex::Real uInf = ddata.reference_velocity & ddata.normal_vec;
        const amrex::Real rho = ddata.density;
        const amrex::Real cT = ddata.thrust_coeff;
        const int npts = ddata.num_force_pts;
        const amrex::Real aeroPressure = 0.5 * uInf * uInf * rho * cT;
        const amrex::Real dr = ddata.dr;
        const auto& normal = ddata.normal_vec;

        for (int ip = 0; ip < npts; ++ip) {
            const amrex::Real r = (ip + 0.5) * dr;
            const amrex::Real rp = r + dr * 0.5;
            const amrex::Real rm = r - dr * 0.5;
            const amrex::Real a = ::amr_wind::utils::pi() * (rp * rp - rm * rm);
            grid.force[ip] = -(aeroPressure * a) * normal;
        }
    }
};

template <>
struct UpdatePosOp<ConstantCt, ActSrcDisk>
{
    void operator()(ConstantCt::DataType&) {}
};

template <>
struct ProcessOutputsOp<ConstantCt, ActSrcDisk>
{
    explicit ProcessOutputsOp<ConstantCt, ActSrcDisk>(ConstantCt::DataType&) {}
    void operator()(ConstantCt::DataType&) {}
    void read_io_options(const utils::ActParser&) {}
    void prepare_outputs(const std::string&) {}
    void write_outputs() {}
};

} // namespace ops
} // namespace actuator
} // namespace amr_wind

#endif /* CONSTANT_CT_OPS_H_ */
