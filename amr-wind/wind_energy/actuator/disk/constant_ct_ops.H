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
        pp.get("num_points", ddata.num_points);
        pp.get("disk_center", ddata.center);
        pp.get("disk_normal", ddata.normal_vec);
        pp.get("epsilon", ddata.epsilon);
        pp.get("diameter", ddata.diameter);

        const amrex::Real radius = ddata.diameter * 0.5;
        ddata.dr = radius / ddata.num_points;

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
    void operator()(ConstantCt::DataType&) {}
};

template <>
struct UpdateVelOp<ConstantCt, ActSrcDisk>
{
    void operator()(ConstantCt::DataType&) {}
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
        const amrex::Real uInf = ddata.reference_velocity;
        const amrex::Real rho = ddata.density;
        const amrex::Real cT = ddata.thrust_coeff;
        const int npts = ddata.num_points;
        const amrex::Real aeroPressure = 0.5 * uInf * uInf * rho * cT;
        const amrex::Real dr = ddata.dr;

        auto* pos = grid.pos.data();
        
        for (int ip = 0; ip < npts; ++ip) {
            const amrex::Real r = (ip+0.5)*dr;
            const amrex::Real rp = r+dr*0.5;
            const amrex::Real rm = r-dr*0.5;
            const amrex::Real a = ::amr_wind::utils::pi()*(rp*rp-rm*rm);
            grid.force(ip) = aeroPressure * a;
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
