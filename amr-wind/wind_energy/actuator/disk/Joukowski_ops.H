#ifndef JOUKOWSKI_OPS_H_
#define JOUKOWSKI_OPS_H_

#include "amr-wind/wind_energy/actuator/disk/Joukowski.H"
#include "amr-wind/wind_energy/actuator/disk/disk_ops.H"
#include "amr-wind/utilities/linear_interpolation.H"

namespace amr_wind {
namespace actuator {
namespace ops {

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

template <>
struct ComputeForceOp<Joukowski, ActSrcDisk>
{
    void operator()(Joukowski::DataType& /*data*/) {}
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
