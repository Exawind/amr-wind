#ifndef FASTIFACE_H
#define FASTIFACE_H

#include "amr-wind/core/ExtSolver.H"
#include "amr-wind/wind_energy/actuator/turbine/fast/fast_wrapper.H"
#include "amr-wind/wind_energy/actuator/turbine/fast/fast_types.H"
#include <map>
#include <vector>

namespace ncutils {
class NCFile;
}

namespace amr_wind {
class CFDSim;
}

namespace exw_fast {

class FastIface : public ::amr_wind::ExtSolver::Register<FastIface>
{
public:
    static std::string identifier() { return "OpenFAST"; }

    explicit FastIface(const ::amr_wind::CFDSim& sim);

    ~FastIface() override;

    void parse_inputs(
        const amr_wind::CFDSim& /*sim*/, const std::string& /*inp_name*/);

    int register_turbine(FastTurbine& data);

    void init_turbine(const int local_id);

    void init_solution(const int local_id);

    void advance_turbine(const int local_id);

    void save_restart(const int local_id);

    int num_local_turbines() const
    {
        return static_cast<int>(m_turbine_data.size());
    }

protected:
    void allocate_fast_turbines();

    void fast_init_turbine(FastTurbine& /*fi*/);

    void fast_restart_turbine(FastTurbine& /*unused*/);

    void fast_replay_turbine(FastTurbine& /*fi*/);

    void prepare_netcdf_file(FastTurbine& /*unused*/);

    void write_velocity_data(const FastTurbine& /*unused*/);

    void read_velocity_data(
        FastTurbine& /*unused*/,
        const ncutils::NCFile& /*unused*/,
        const size_t tid);

    //! Global to local index lookup map
    std::map<int, int> m_turbine_map;

    std::vector<FastTurbine*> m_turbine_data;

    std::string m_output_dir{"fast_velocity_data"};

    double m_dt_cfd{0.0};

    // Not using amrex::Real. Instead use explicit types from FAST API
    double m_start_time{0.0};
    double m_stop_time{-1.0};
    SimMode m_sim_mode{SimMode::init};

    int m_num_sc_inputs{0};
    int m_num_sc_outputs{0};

#ifdef AMR_WIND_FAST_USE_SCDX
    int m_num_sc_inputs_glob{0};
    float m_init_sc_inputs_glob{0.0};
    float m_init_sc_inputs_turbine{0.0};
#endif

    bool m_is_initialized{false};
};

} // namespace exw_fast

#endif /* FASTIFACE_H */
