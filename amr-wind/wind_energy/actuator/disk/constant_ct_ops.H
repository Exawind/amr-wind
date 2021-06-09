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
        pp.get("disk_normal", ddata.normal);
        pp.get("epsilon", ddata.epsilon);
        pp.get("diameter", ddata.diameter);
        // TODO think about center of disk
        const amrex::Real radius = ddata.diameter * 0.5;
        ddata.dr = radius / ddata.num_points;

        // ensure normal is normalized
        ddata.normal = ddata.normal.normalize();
        const vs::Vector& n = ddata.normal;
        // get a vector in the plane of the disk
        vs::Vector planeVec = vs::Vector::one();
        const auto unit = vs::Tensor::I();

        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            std::abs((unit.z() & n) - vs::mag_sqr(n)) >
                vs::DTraits<amrex::Real>::eps(),
            "A constant Ct disk normal is too close to vertical. Vertical axis "
            "wind turbines are not supported with this model.");

        ::amr_wind::utils::cross_prod(
            n.data(), unit.z().data(), planeVec.data());

        // clang-format off
        const amrex::Real nl = ddata.epsilon * 3.0; // length scale in normal dir
        const amrex::Real dl = radius + ddata.dr*2.0; // length scale in plane of disk

        const auto& cc = ddata.center;
        const amrex::Real xLen = (unit.x()&n)*nl + (unit.x()&planeVec)*dl;
        const amrex::Real yLen = (unit.y()&n)*nl + (unit.y()&planeVec)*dl;
        const amrex::Real zLen = (unit.z()&n)*nl + (unit.z()&planeVec)*dl;

        info.bound_box = amrex::RealBox(
          cc.x()-xLen, cc.y()-yLen, cc.z()-zLen, 
          cc.x()+xLen, cc.y()+yLen, cc.z()+zLen
        );
        // clang-format on
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

// decide if we want force defined on the grid.  could
// make grid just for velocity probing and then have
// seperate points for the force that live in DataType.meta
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
        for (int ip = 0; ip < npts; ++ip) {
            // get radius at point
            // get ring area
            const amrex::Real area = 1.0;
            grid.force(ip) = aeroPressure * area;
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
