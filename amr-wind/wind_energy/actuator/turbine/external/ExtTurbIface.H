#ifndef EXTTURBIFACE_H
#define EXTTURBIFACE_H

#include "amr-wind/core/ExtSolver.H"
#include "amr-wind/wind_energy/actuator/turbine/external/external_base_types.H"
#include <map>
#include <vector>

namespace ncutils {
class NCFile;
}

namespace amr_wind {
class CFDSim;
}

namespace exw_fast {

template <typename SolverTurbine>
class ExtTurbIface
    : public ::amr_wind::ExtSolver::Register<ExtTurbIface<SolverTurbine>>
{
public:
    static std::string identifier() { return "OpenFAST"; }

    explicit ExtTurbIface(const ::amr_wind::CFDSim& sim);

    ~ExtTurbIface() override;

    void parse_inputs(
        const amr_wind::CFDSim& /*sim*/, const std::string& /*inp_name*/);

    int register_turbine(SolverTurbine& data);

    void init_turbine(const int local_id);

    void init_solution(const int local_id);

    void get_hub_stats(const int local_id);

    void advance_turbine(const int local_id);

    void save_restart(const int local_id);

    int num_local_turbines() const
    {
        return static_cast<int>(m_turbine_data.size());
    }

protected:
    void allocate_fast_turbines();

    void fast_init_turbine(SolverTurbine& /*fi*/);

    void fast_restart_turbine(SolverTurbine& /*fi*/);

    void fast_replay_turbine(SolverTurbine& /*fi*/);

    void prepare_netcdf_file(SolverTurbine& /*unused*/);

    void write_velocity_data(const SolverTurbine& /*unused*/);

    static void read_velocity_data(
        SolverTurbine& /*unused*/,
        const ncutils::NCFile& /*unused*/,
        const size_t tid);

    //! Global to local index lookup map
    std::map<int, int> m_turbine_map;

    std::vector<SolverTurbine*> m_turbine_data;

    std::string m_output_dir{"fast_velocity_data"};

    double m_dt_cfd{0.0};

    // Not using amrex::Real. Instead use explicit types from FAST API
    double m_start_time{0.0};
    double m_stop_time{-1.0};
    ::exw_fast::SimMode m_sim_mode{::exw_fast::SimMode::init};

    int m_num_sc_inputs{0};
    int m_num_sc_outputs{0};

#if OPENFAST_VERSION_MAJOR == 4
    int m_inflow_type{2};
#endif

    int m_num_sc_inputs_glob{0};
    float m_init_sc_inputs_glob{0.0};
    float m_init_sc_inputs_turbine{0.0};

    bool m_is_initialized{false};
};

// General implementations of some functions
template <typename SolverTurbine>
ExtTurbIface<SolverTurbine>::ExtTurbIface(const amr_wind::CFDSim& /*unused*/)
{}

} // namespace exw_fast

#endif /* EXTTURBIFACE_H */
