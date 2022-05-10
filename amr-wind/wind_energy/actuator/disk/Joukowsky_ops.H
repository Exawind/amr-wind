#ifndef Joukowsky_OPS_H_
#define Joukowsky_OPS_H_

#include "amr-wind/wind_energy/actuator/disk/Joukowsky.H"
#include "amr-wind/wind_energy/actuator/disk/disk_ops.H"
#include "amr-wind/utilities/linear_interpolation.H"

namespace amr_wind {
namespace actuator {
namespace ops {

namespace joukowsky {
void parse_and_gather_params(const utils::ActParser& pp, JoukowskyData& data);
inline void
set_current_angular_velocity(JoukowskyData& data, const amrex::Real uInfSqr)
{
    const amrex::Real uInfMag = std::sqrt(uInfSqr);
    data.current_angular_velocity = ::amr_wind::interp::linear(
        data.table_velocity, data.angular_velocity, uInfMag);
}
} // namespace joukowsky

template <>
struct ReadInputsOp<Joukowsky, ActSrcDisk>
{
    void operator()(Joukowsky::DataType& data, const utils::ActParser& pp)
    {
        auto& meta = data.meta();
        joukowsky::parse_and_gather_params(pp, meta);
        base::do_parse_based_computations<Joukowsky>(data);
    }
};

template <>
struct InitDataOp<Joukowsky, ActSrcDisk>
{
    void operator()(typename Joukowsky::DataType& data)
    {
        auto& grid = data.grid();
        auto& meta = data.meta();

        // only resize the members we are going to use
        grid.pos.resize(meta.num_force_pts);
        grid.force.resize(meta.num_force_pts);
        grid.vel.resize(meta.num_vel_pts);
        grid.vel_pos.resize(meta.num_vel_pts);
        meta.tip_correction.resize(meta.num_vel_pts_r);
        meta.root_correction.resize(meta.num_vel_pts_r);
        std::fill(meta.tip_correction.begin(), meta.tip_correction.end(), 1.0);
        std::fill(
            meta.root_correction.begin(), meta.root_correction.end(), 1.0);

        const auto& sVec = meta.sample_vec;
        const auto& nVec = meta.normal_vec;

        // force points
        base::compute_disk_points(meta, grid.pos, nVec, 0, 0);
        // velocity points upstream
        base::compute_disk_points(
            meta, grid.vel_pos, sVec, 0, meta.diameters_to_sample);
        // velocity points at the disk
        base::compute_disk_points(
            meta, grid.vel_pos, nVec, meta.num_vel_pts / 2, 0);
    }
};

/**
  @brief Compute Forces following procedures of Sorenson 2020

  Follow the procedure for computing body forces and applying tip/root
  corrections using a Joukowsky disk as outlined in:

  Sørensen, Jens Nørkær, et al. "Analytical body forces in numerical actuator
  disc model of wind turbines." Renewable Energy 147 (2020): 2259-2271.
 */
template <>
struct ComputeForceOp<Joukowsky, ActSrcDisk>
{
    void operator()(Joukowsky::DataType& data)
    {
        // Equation comments refer to Sorenson 2020 (full reference above)
        auto& grid = data.grid();
        auto& ddata = data.meta();
        const amrex::Real machine_eps = vs::DTraits<amrex::Real>::eps();
        // step 1: compute Ct, a1 and a2 coefficients (eq 16)
        const amrex::Real uInfSqr = base::compute_reference_velocity_sqr(ddata);
        const amrex::Real U_ref = std::max(std::sqrt(uInfSqr), machine_eps);
        base::set_thrust_coefficient(ddata, uInfSqr);
        joukowsky::set_current_angular_velocity(ddata, uInfSqr);

        const amrex::Real Ct = ddata.current_ct;

        amrex::Real a1 = 0.0;
        amrex::Real a2 = 0.0;

        for (int ip = 0; ip < ddata.num_vel_pts_r; ip++) {
            const amrex::Real r = (ip + 0.5) * ddata.dr;

            a1 += std::pow(ddata.tip_correction[ip], 2.0) *
                  std::pow(ddata.root_correction[ip], 2.0) / r * ddata.dr;

            a2 += ddata.tip_correction[ip] * ddata.root_correction[ip] * r *
                  ddata.dr;
        }

        // step 2: determine the circulation (q0) from a1, a2 and Ct (eq 17)

        const amrex::Real lambda =
            0.5 * ddata.diameter * ddata.current_angular_velocity / U_ref;

        const amrex::Real q0 =
            (std::sqrt(16.0 * lambda * lambda * a2 * a2 + 8.0 * a1 * Ct) -
             4.0 * lambda * a2) /
            std::max(4.0 * a1, machine_eps);

        // step 3: compute normal force (fz) and azimuthal force (f_theta) (eq
        // 13)

        VecSlice disk_velocity = ::amr_wind::utils::slice(
            grid.vel, ddata.num_force_pts, ddata.num_force_pts);
        int ip = 0;
        for (int i = 0; i < ddata.num_vel_pts_r; i++) {
            const amrex::Real& F = ddata.tip_correction[i];
            const amrex::Real& g = ddata.root_correction[i];

            const amrex::Real x = std::max(
                (i + 0.5) * ddata.dr / (0.5 * ddata.diameter), machine_eps);

            const amrex::Real f_z =
                q0 * g * F / x * (lambda * x + 0.5 * q0 * g * F / x);

            const auto point_on_coplanar_vector = grid.pos[i];

            for (int j = 0; j < ddata.num_vel_pts_t; j++, ip++) {
                const auto point_current = grid.pos[ip];

                // compute angle between points on the coplanar vector and the
                // current points
                const auto distance = utils::delta_pnts_cyl(
                    ddata.center, ddata.normal_vec, point_on_coplanar_vector,
                    point_current);

                // TODO add sign convention for rotation (+/- RPM)
                const auto angle = ::amr_wind::utils::degrees(-distance[1]);

                const auto theta_vec = vs::quaternion(ddata.normal_vec, angle) &
                                       ddata.coplanar_vec;

                const amrex::Real u_disk_ij =
                    ddata.normal_vec & disk_velocity[ip];

                const amrex::Real f_theta = u_disk_ij / U_ref * q0 * g * F / x;

                grid.force[ip] =
                    (f_z * ddata.normal_vec + f_theta * theta_vec) *
                    ddata.density * uInfSqr;
            }
        }

        // step 5: compute moment and power (eqs 18 and 19)
        // step 6: compute cp
    }
};

template <>
struct ProcessOutputsOp<Joukowsky, ActSrcDisk>
{
private:
    Joukowsky::DataType& m_data;
    //! Path to the output directory (specified by Actuator physics class)
    std::string m_out_dir{""};

    //! NetCDF output filename for this turbine
    std::string m_nc_filename{""};

    //! Output frequency (specified in input file)
    int m_out_freq{10};

public:
    // cppcheck-suppress constParameter
    explicit ProcessOutputsOp<Joukowsky, ActSrcDisk>(Joukowsky::DataType& data)
        : m_data(data)
    {}
    void operator()(Joukowsky::DataType& /*data*/) {}
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

#endif /* Joukowsky_OPS_H_ */
