#ifndef CONSTANT_CT_OPS_H_
#define CONSTANT_CT_OPS_H_

#include "amr-wind/wind_energy/actuator/disk/ConstantCt.H"
#include "amr-wind/wind_energy/actuator/actuator_ops.H"
#include "amr-wind/utilities/tensor_ops.H"

namespace amr_wind {
namespace actuator {
namespace ops {

template <>
struct ReadInputsOp<ConstantCt, ActSrcDisk>
{
    void operator()(ConstantCt::DataType& data, const utils::ActParser& pp)
    {
        auto& ddata = data.meta();
        auto& info = data.info();
        pp.get("num_force_points", ddata.num_force_pts);
        // eventually we will change velocity sampling location
        // pp.get("num_vel_points_r", ddata.num_vel_points_r);
        // pp.get("num_vel_points_t", ddata.num_vel_points_t);
        // pp.get("diameters_to_sample_velocity", ddata.diameters_to_sample);
        pp.get("disk_center", ddata.center);
        pp.get("disk_normal", ddata.normal_vec);
        pp.get("epsilon", ddata.epsilon);
        pp.get("diameter", ddata.diameter);
        pp.get("thrust_coeff", ddata.thrust_coeff);
        pp.query("density", ddata.density);

        // ddata.num_vel_pts = ddata.num_vel_pts_r * ddata.num_vel_pts_t;
        ddata.num_vel_pts = ddata.num_force_pts;
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            ddata.num_vel_pts > 0,
            "num_vel_points_r and num_vel_points_t must both be >=1");

        const amrex::Real radius = ddata.diameter * 0.5;
        ddata.dr = radius / ddata.num_force_pts;

        // ensure normal is normalized
        ddata.normal_vec.normalize();
        const vs::Vector& norm = ddata.normal_vec;

        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            vs::mag_sqr(vs::Vector::khat()^norm) <
                vs::DTraits<amrex::Real>::eps(),
            "A constant Ct disk normal is too close to vertical.");
        
        // compute a coplanar vector that resides in the same plane as the disk
        // we will use this vector for the bounding box and the force point locations
        ddata.coplanar_vec = norm ^ vs::Vector::khat();
        ddata.coplanar_vec.normalize();
        const vs::Vector& cVec = ddata.coplanar_vec;

        const amrex::Real nl = ddata.epsilon * 3.0; // length scale in normal dir
        const amrex::Real dl = radius + ddata.dr*2.0; // length scale in plane of disk

        const auto& cc = ddata.center;
        const auto dvec = norm * nl + cVec * dl + vs::Vector::khat() * dl;
        const auto p1 = cc - dvec;  // front
        const auto p2 = cc + dvec;  // back
        info.bound_box = amrex::RealBox(
            amrex::min(p1.x(), p2.x()),
            amrex::min(p1.y(), p2.y()),
            amrex::min(p1.z(), p2.z()),
            amrex::max(p1.x(), p2.x()),
            amrex::max(p1.y(), p2.y()),
            amrex::max(p1.z(), p2.z())
        );

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
        // const auto& nVec = meta.normal_vec;
        const auto& cc = meta.center;
        {
            const auto& dr = meta.dr;
            for (int i = 0; i < meta.num_force_pts; ++i) {
                grid.pos[i] = cc + (i + 0.5) * dr * cVec;
                grid.vel_pos[i] = grid.pos[i];
            }
        }
        /* TODO compute a disk x diameters away that can be averaged to get
        a reference velocity and density
        {
            const int nR = meta.num_vel_pts_r;
            const int nT = meta.num_vel_pts_t;
            const amrex::Real dr = meta.diameter / nR;
            const amrex::Real dt = ::amr_wind::utils::two_pi() / nT;
            const amrex::Real dU = meta.diameters_to_sample;
            int k = 0;
            for (int i = 0; i < nR; i++) {
                const amrex::Real r = dr * (i + 0.5);
                for (int j = 0; j < nT; j++, k++) {
                    const amrex::Real theta = j * dt;
                    const vs::Vector relativePosition = {
                        r * std::cos(theta), r * std::sin(theta), dU};
                    grid.vel_pos[k] = cc + dU*nVec +
                }
            }
        }
        */
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

        auto* pos = grid.pos.data();
        
        for (int ip = 0; ip < npts; ++ip) {
            const amrex::Real r = (ip+0.5)*dr;
            const amrex::Real rp = r+dr*0.5;
            const amrex::Real rm = r-dr*0.5;
            const amrex::Real a = ::amr_wind::utils::pi()*(rp*rp-rm*rm);
            grid.force[ip] = (aeroPressure * a) * normal;
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
