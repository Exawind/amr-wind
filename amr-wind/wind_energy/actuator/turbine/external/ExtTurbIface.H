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

    void do_turbine_step(int& tid);

    void write_turbine_checkpoint(int& tid);

    void save_restart(const int local_id);

    int num_local_turbines() const
    {
        return static_cast<int>(m_turbine_data.size());
    }

protected:
    void allocate_ext_turbines();

    void ext_init_turbine(SolverTurbine& /*fi*/);

    void ext_restart_turbine(SolverTurbine& /*fi*/);

    void ext_replay_turbine(SolverTurbine& /*fi*/);

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

template <typename SolverTurbine>
void ExtTurbIface<SolverTurbine>::advance_turbine(const int local_id)
{
    BL_PROFILE("amr-wind::FastIface::advance_turbine");
    AMREX_ASSERT(local_id < static_cast<int>(m_turbine_data.size()));

    auto& fi = *m_turbine_data[local_id];
    AMREX_ASSERT(!fi.is_solution0);
    {
        const auto& tmax = fi.stop_time;
        const auto& telapsed = (fi.time_index + fi.num_substeps) * fi.dt_fast;
        if (telapsed > (tmax + 1.0e-8)) {
            // clang-format off
            amrex::OutStream()
                << "\nWARNING: FastIface:\n"
                << "  Elapsed simulation time will exceed max "
                << "time set for OpenFAST"
                << std::endl << std::endl;
            // clang-format on
        }
    }

    write_velocity_data(fi);
    for (int i = 0; i < fi.num_substeps; ++i, ++fi.time_index) {
        do_turbine_step(fi.tid_local);
    }

    if (fi.chkpt_interval > 0 &&
        (fi.time_index / fi.num_substeps) % fi.chkpt_interval == 0) {
        write_turbine_checkpoint(fi.tid_local);
    }
}

template <typename SolverTurbine>
void ExtTurbIface<SolverTurbine>::init_turbine(const int local_id)
{
    AMREX_ALWAYS_ASSERT(local_id < static_cast<int>(m_turbine_data.size()));
    if (!m_is_initialized) {
        allocate_ext_turbines();
    }
    auto& fi = *m_turbine_data[local_id];

    switch (fi.sim_mode) {
    case ::exw_fast::SimMode::init: {
        ext_init_turbine(fi);
        prepare_netcdf_file(fi);
        break;
    }

    case ::exw_fast::SimMode::replay: {
        ext_init_turbine(fi);
        ext_replay_turbine(fi);
        break;
    }

    case ::exw_fast::SimMode::restart: {
        ext_restart_turbine(fi);
        prepare_netcdf_file(fi);
        break;
    }
    }
}

} // namespace exw_fast

#endif /* EXTTURBIFACE_H */
