#ifndef UNIFORM_CT_OPS_H_
#define UNIFORM_CT_OPS_H_

#include "amr-wind/wind_energy/actuator/disk/UniformCt.H"
#include "amr-wind/wind_energy/actuator/actuator_ops.H"
#include "amr-wind/utilities/tensor_ops.H"
#include "amr-wind/utilities/linear_interpolation.H"

namespace amr_wind {
namespace actuator {
namespace disk {
// refactor these when we add more disk types later
void prepare_netcdf_file(
    const std::string&, const UniformCtData&, const ActInfo&, const ActGrid&);

void write_netcdf(
    const std::string&,
    const UniformCtData&,
    const ActInfo&,
    const ActGrid&,
    const amrex::Real);
} // namespace disk

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

void required_parameters(UniformCt::MetaType& meta, const utils::ActParser& pp);

void optional_parameters(UniformCt::MetaType& meta, const utils::ActParser& pp);

void check_for_parse_conflicts(const utils::ActParser& pp);

void compute_and_normalize_coplanar_vector(UniformCt::MetaType& meta);

void final_checks(const UniformCt::MetaType& meta);

amrex::RealBox compute_bounding_box(const UniformCt::MetaType& meta);

void do_parse_based_computations(UniformCt::DataType& data);

void compute_disk_points(
    UniformCt::MetaType& meta,
    VecList& points,
    const vs::Vector& cylAxis,
    const int offset,
    const double dOffset);

template <>
struct ReadInputsOp<UniformCt, ActSrcDisk>
{
    void operator()(UniformCt::DataType& data, const utils::ActParser& pp)
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
struct InitDataOp<UniformCt, ActSrcDisk>
{
    void operator()(UniformCt::DataType& data)
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
        compute_disk_points(
            meta, grid.vel_pos, sVec, 0, meta.diameters_to_sample);
        compute_disk_points(meta, grid.vel_pos, nVec, meta.num_vel_pts / 2, 0);
    }
};

template <>
struct UpdateVelOp<UniformCt, ActSrcDisk>
{
    void operator()(UniformCt::DataType& data)
    {
        auto& meta = data.meta();
        const auto& grid = data.grid();
        meta.reference_velocity = {0.0, 0.0, 0.0};
        meta.disk_velocity = {0.0, 0.0, 0.0};
        auto& refVel = meta.reference_velocity;
        auto& diskVel = meta.disk_velocity;
        const int np = meta.num_vel_pts / 2;
        for (int i = 0; i < np; i++) {
            refVel = refVel + grid.vel[i];
            diskVel = diskVel + grid.vel[i + np];
        }
        refVel /= np;
        diskVel /= np;
    }
};

// compute the total force over a given radial section
// this will then be spread uniformly over that area in the
// projection step
template <>
struct ComputeForceOp<UniformCt, ActSrcDisk>
{
    void operator()(UniformCt::DataType& data)
    {
        auto& grid = data.grid();
        auto& ddata = data.meta();
        const auto& normal = ddata.normal_vec;
        const amrex::Real rho = ddata.density;
        // only use veloocity compnents normal to disk
        // square the velocity by componet first so normal direction signs will
        // cancel each other
        const amrex::Real uInfSqr =
            (ddata.reference_velocity * ddata.reference_velocity) & normal;

        ddata.current_ct = ::amr_wind::interp::linear(
            ddata.table_velocity, ddata.thrust_coeff,
            std::sqrt(std::abs(uInfSqr)));

        const int npts = ddata.num_force_pts;
        const amrex::Real aeroPressure = 0.5 * uInfSqr * rho * ddata.current_ct;
        const amrex::Real dr = ddata.dr;

        for (int ip = 0; ip < npts; ++ip) {
            const amrex::Real r = (ip + 0.5) * dr;
            const amrex::Real rp = r + dr * 0.5;
            const amrex::Real rm = r - dr * 0.5;
            const amrex::Real a = ::amr_wind::utils::pi() * (rp * rp - rm * rm);
            // canceling of signes ensures this term will always be negative
            grid.force[ip] = -(aeroPressure * a) * normal;
        }
    }
};

template <>
struct UpdatePosOp<UniformCt, ActSrcDisk>
{
    void operator()(UniformCt::DataType&) {}
};

template <>
struct ProcessOutputsOp<UniformCt, ActSrcDisk>
{
private:
    UniformCt::DataType& m_data;
    //! Path to the output directory (specified by Actuator physics class)
    std::string m_out_dir{""};

    //! NetCDF output filename for this turbine
    std::string m_nc_filename{""};

    //! Output frequency (specified in input file)
    int m_out_freq{10};

public:
    explicit ProcessOutputsOp<UniformCt, ActSrcDisk>(UniformCt::DataType& data)
        : m_data(data)
    {}
    void operator()(UniformCt::DataType&) {}
    void read_io_options(const utils::ActParser& pp)
    {
        pp.query("output_frequency", m_out_freq);
    }
    void prepare_outputs(const std::string& out_dir)
    {
        m_nc_filename = out_dir + "/" + m_data.info().label + ".nc";
        disk::prepare_netcdf_file(
            m_nc_filename, m_data.meta(), m_data.info(), m_data.grid());
    }
    void write_outputs()
    {
        const auto& time = m_data.sim().time();
        const int tidx = time.time_index();
        if (tidx % m_out_freq != 0) return;

        disk::write_netcdf(
            m_nc_filename, m_data.meta(), m_data.info(), m_data.grid(),
            time.new_time());
    }
};

} // namespace ops
} // namespace actuator
} // namespace amr_wind

#endif /* UNIFORM_CT_OPS_H_ */
