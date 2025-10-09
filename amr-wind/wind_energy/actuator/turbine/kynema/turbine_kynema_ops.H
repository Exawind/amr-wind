#ifndef TURBINE_KYNEMA_OPS_H
#define TURBINE_KYNEMA_OPS_H

#include "amr-wind/wind_energy/actuator/turbine/kynema/TurbineKynema.H"
#include "amr-wind/wind_energy/actuator/turbine/external/turbine_external_ops.H"

namespace amr_wind::actuator::ops {

template <typename SrcTrait>
struct ReadInputsOp<TurbineKynema, SrcTrait>
{
    void operator()(TurbineKynema::DataType& data, const utils::ActParser& pp)
    {
        read_ops<TurbineKynema::DataType>(data, pp);

        auto& tf = data.meta().ext_data;

        // Only one available for now
        const std::string sim_mode{"init"};

        if (sim_mode == "init") {
            tf.sim_mode = ::ext_turb::SimMode::init;
            amrex::Print() << "Initializing turbine:" << tf.tlabel << std::endl;
        } else if (sim_mode == "replay") {
            tf.sim_mode = ::ext_turb::SimMode::replay;
            amrex::Print() << "Replaying turbine:" << tf.tlabel << std::endl;
        } else if (sim_mode == "restart") {
            tf.sim_mode = ::ext_turb::SimMode::restart;
            amrex::Print() << "Restarting turbine:" << tf.tlabel << std::endl;
        } else {
            amrex::Abort(
                "Actuator: Invalid Kynema simulation mode: " + sim_mode);
        }

        if (tf.sim_mode == ::ext_turb::SimMode::restart) {
            pp.get("kynema_restart_file", tf.checkpoint_file);
        } else {
            pp.get("kynema_input_file", tf.input_file);
        }
    }
};

template <>
inline void determine_influenced_procs<TurbineKynema>(
    typename TurbineKynema::DataType& data)
{
    external_determine_influenced_procs<TurbineKynema::DataType>(data);
}

template <>
inline void determine_root_proc<TurbineKynema>(
    typename TurbineKynema::DataType& data, amrex::Vector<int>& act_proc_count)
{
    external_determine_root_proc<
        TurbineKynema::DataType, ::ext_turb::KynemaTurbine,
        ::ext_turb::KynemaSolverData>(data, act_proc_count);
}

template <typename SrcTrait>
struct InitDataOp<TurbineKynema, SrcTrait>
{
    void operator()(TurbineKynema::DataType& data)
    {
        BL_PROFILE("amr-wind::InitDataOp<TurbineKynema>");
        init_data_op<TurbineKynema::DataType>(data);
    }
};

template <typename SrcTrait>
struct UpdatePosOp<TurbineKynema, SrcTrait>
{
    void operator()(typename TurbineKynema::DataType& data)
    {
        BL_PROFILE("amr-wind::actuator::UpdatePosOp<TurbineKynema>");
        update_pos_op<TurbineKynema::DataType>(data);
    }
};

template <typename SrcTrait>
struct UpdateVelOp<TurbineKynema, SrcTrait>
{
    void operator()(typename TurbineKynema::DataType& data)
    {
        BL_PROFILE("amr-wind::actuator::UpdateVelOp<TurbineKynema>");
        update_vel_op<TurbineKynema::DataType>(data);
    }
};

template <typename SrcTrait>
struct ComputeForceOp<TurbineKynema, SrcTrait>
{
    void operator()(typename TurbineKynema::DataType& data)
    {
        BL_PROFILE("amr-wind::actuator::ComputeForceOp<TurbineKynema>");
        compute_force_op<TurbineKynema::DataType>(data);
    }
};

template <typename SrcTrait>
struct ProcessOutputsOp<TurbineKynema, SrcTrait>
{
private:
    typename TurbineKynema::DataType& m_data;

    //! Path to the output directory (specified by Actuator physics class)
    std::string m_out_dir;

    //! NetCDF output filename for this turbine
    std::string m_nc_filename;

    //! Output frequency (specified in input file)
    int m_out_freq{10};

public:
    explicit ProcessOutputsOp(typename TurbineKynema::DataType& data)
        : m_data(data)
    {}

    void read_io_options(const utils::ActParser& pp)
    {
        pp.query("output_frequency", m_out_freq);
    }

    void prepare_outputs(const std::string& out_dir)
    {
        m_nc_filename = out_dir + "/" + m_data.info().label + ".nc";
        utils::prepare_netcdf_file(
            m_nc_filename, m_data.meta(), m_data.info(), m_data.grid());
    }

    void write_outputs()
    {
        const auto& time = m_data.sim().time();
        const int tidx = time.time_index();
        if ((m_out_freq > 0) && (tidx % m_out_freq != 0)) {
            return;
        }

        utils::write_netcdf(
            m_nc_filename, m_data.meta(), m_data.info(), m_data.grid(),
            time.new_time());
    }
};

} // namespace amr_wind::actuator::ops

#endif /* TURBINE_KYNEMA_OPS_H */
