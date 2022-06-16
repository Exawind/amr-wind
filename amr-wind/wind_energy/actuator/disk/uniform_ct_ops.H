#ifndef UNIFORM_CT_OPS_H_
#define UNIFORM_CT_OPS_H_

#include "amr-wind/wind_energy/actuator/disk/disk_ops.H"
#include "amr-wind/wind_energy/actuator/disk/UniformCt.H"
#include "amr-wind/utilities/linear_interpolation.H"
namespace amr_wind {
namespace actuator {
namespace ops {
namespace uniformct {

void parse_and_gather_params(const utils::ActParser& pp, UniformCtData& data);

} // namespace uniformct
template <>
struct ReadInputsOp<UniformCt, ActSrcDisk>
{
    void operator()(UniformCt::DataType& data, const utils::ActParser& pp)
    {
        // TODO get density at probe location
        auto& meta = data.meta();
        uniformct::parse_and_gather_params(pp, meta);
        base::do_parse_based_computations<UniformCt>(data);
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
        const amrex::Real uInfSqr = base::compute_reference_velocity_sqr(ddata);
        base::set_thrust_coefficient(ddata, uInfSqr);

        const int npts = ddata.num_force_pts;
        const amrex::Real rho = ddata.density;
        const amrex::Real aeroPressure = 0.5 * uInfSqr * rho * ddata.current_ct;
        const amrex::Real dr = ddata.dr;

        for (int ip = 0; ip < npts; ++ip) {
            const amrex::Real r = (ip + 0.5) * dr;
            const amrex::Real rp = r + dr * 0.5;
            const amrex::Real rm = r - dr * 0.5;
            const amrex::Real a = ::amr_wind::utils::pi() * (rp * rp - rm * rm);
            // disk force should always point in direction of the normal
            grid.force[ip] = (aeroPressure * a) * normal;
        }
    }
};

template <>
struct InitDataOp<UniformCt, ActSrcDisk>
{
    void operator()(typename UniformCt::DataType& data)
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
        base::compute_disk_points(
            meta, grid.vel_pos, sVec, 0, meta.diameters_to_sample);
        base::compute_disk_points(
            meta, grid.vel_pos, nVec, meta.num_vel_pts / 2, 0);
    }
};
template <>
struct ProcessOutputsOp<UniformCt, ActSrcDisk>
{
private:
    // cppcheck-suppress uninitMemberVarPrivate
    UniformCt::DataType& m_data;
    //! Path to the output directory (specified by Actuator physics class)
    std::string m_out_dir;

    //! NetCDF output filename for this turbine
    std::string m_nc_filename;

    //! Output frequency (specified in input file)
    int m_out_freq{10};

public:
    // cppcheck-suppress constParameter
    explicit ProcessOutputsOp<UniformCt, ActSrcDisk>(UniformCt::DataType& data)
        : m_data(data)
    {}
    void operator()(UniformCt::DataType& /*unused*/) {}
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
        if (tidx % m_out_freq != 0) {
            return;
        }

        disk::write_netcdf(
            m_nc_filename, m_data.meta(), m_data.info(), m_data.grid(),
            time.new_time());
    }
};
} // namespace ops
} // namespace actuator
} // namespace amr_wind

#endif /* UNIFORM_CT_OPS_H_ */
