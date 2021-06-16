#ifndef CONSTANT_CT_OPS_H_
#define CONSTANT_CT_OPS_H_

#include "amr-wind/wind_energy/actuator/disk/ConstantCt.H"
#include "amr-wind/wind_energy/actuator/actuator_ops.H"
#include "amr-wind/utilities/tensor_ops.H"

namespace amr_wind {
namespace actuator {
namespace ops {

void collect_parse_conflicts(
    const utils::ActParser& pp,
    const std::string& p1,
    const std::string& p2,
    std::ostringstream& ss);
void collect_parse_dependencies(
    const utils::ActParser& pp,
    const std::string& p1,
    const std::string& p2,
    std::ostringstream& ss);

void required_parameters(
    ConstantCt::MetaType& meta, const utils::ActParser& pp);

void optional_parameters(
    ConstantCt::MetaType& meta, const utils::ActParser& pp);

void check_for_parse_conflicts(const utils::ActParser& pp);

void compute_and_normalize_vectors(ConstantCt::MetaType& meta);

void final_checks(const ConstantCt::MetaType& meta);

amrex::RealBox compute_bounding_box(const ConstantCt::MetaType& meta);

void do_parse_based_computations(ConstantCt::DataType& data);

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
