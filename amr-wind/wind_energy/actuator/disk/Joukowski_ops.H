#ifndef JOUKOWSKI_OPS_H_
#define JOUKOWSKI_OPS_H_

#include "amr-wind/wind_energy/actuator/disk/Joukowski.H"
#include "amr-wind/wind_energy/actuator/disk/disk_ops.H"
#include "amr-wind/utilities/linear_interpolation.H"

namespace amr_wind {
namespace actuator {
namespace ops {

namespace joukowsi {}

template <>
struct ReadInputsOp<Joukowski, ActSrcDisk>
{
    void operator()(Joukowski::DataType& data, const utils::ActParser& pp)
    {
        auto& meta = data.meta();
        base::parse_and_gather_params<Joukowski::MetaType>(pp, meta);
        base::final_checks(meta);
        base::do_parse_based_computations<Joukowski>(data);
        amrex::Abort("JoukoskiDisk is under development and not functional.");
    }
};

/**
  @brief Compute Forces following procedures of Sorenson 2020

  Follow the procedure for computing body forces and applying tip/root
  corrections using a Joukowski disk as outlined in:

  Sørensen, Jens Nørkær, et al. "Analytical body forces in numerical actuator
  disc model of wind turbines." Renewable Energy 147 (2020): 2259-2271.
 */
template <>
struct ComputeForceOp<Joukowski, ActSrcDisk>
{
    void operator()(Joukowski::DataType& data)
    {
        // Equation comments refer to Sorenson 2020 (full reference above)
        auto& grid = data.grid();
        auto& ddata = data.meta();
        // step 1: compute Ct, a1 and a2 coefficients (eq 16)
        const amrex::Real uInfSqr = base::compute_reference_velocity_sqr(ddata);
        base::set_thrust_coefficient(ddata, uInfSqr);
        const amrex::Real Ct = ddata.current_ct;

        amrex::Real a1 = 0.0;
        amrex::Real a2 = 0.0;

        const int npts = ddata.num_force_pts;
        const amrex::Real dx = 2.0 * ddata.dr / ddata.diameter;
        for (int ip = 0; ip < npts; ip++) {
            const amrex::Real x = (ip + 0.5) * dx;
            a1 += std::pow(ddata.tip_correction[ip], 2.0) *
                  std::pow(ddata.root_correction[ip], 2.0) / x * dx;
            a2 += ddata.tip_correction[ip] * ddata.root_correction[ip] * x * dx;
        }

        // step 2: determine the circulation (q0) from a1, a2 and Ct (eq 17)
        const amrex::Real lambda = ddata.tip_speed_ratio;
        const amrex::Real q0 =
            (std::sqrt(16.0 * lambda * lambda * a2 * a2 + 8.0 * a1 * Ct) -
             4.0 * lambda * a2) /
            std::max(4.0 * a1, 1e-12);

        // step 3: compute normal force (fz) and azimuthal force (f_theta) (eq
        // 13)
        for (int ip = 0; ip < npts; ip++) {
            const amrex::Real x = (ip + 0.5) * dx;
            const amrex::Real& F = ddata.tip_correction[ip];
            const amrex::Real& g = ddata.root_correction[ip];
            // loop over azimuthal points TODO(psakiev) make sure
            // num_force_pts_theta can be greater than 1 for LinearBasis
            // spreading
        }

        // step 5: compute moment and power (eqs 18 and 19)
        // step 6: compute cp
    }
};

template <>
struct ProcessOutputsOp<Joukowski, ActSrcDisk>
{
private:
    Joukowski::DataType& m_data;
    //! Path to the output directory (specified by Actuator physics class)
    std::string m_out_dir{""};

    //! NetCDF output filename for this turbine
    std::string m_nc_filename{""};

    //! Output frequency (specified in input file)
    int m_out_freq{10};

public:
    // cppcheck-suppress constParameter
    explicit ProcessOutputsOp<Joukowski, ActSrcDisk>(Joukowski::DataType& data)
        : m_data(data)
    {}
    void operator()(Joukowski::DataType&) {}
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

#endif /* JOUKOWSKI_OPS_H_ */
