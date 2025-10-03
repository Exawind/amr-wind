#ifndef TURBINEEXTERNAL_H
#define TURBINEEXTERNAL_H

#include "amr-wind/wind_energy/actuator/turbine/turbine_types.H"
#include "amr-wind/wind_energy/actuator/turbine/external/ExtTurbIface.H"
#include "amr-wind/core/ExtSolver.H"

namespace amr_wind::actuator {

template <typename SolverTurbine, typename SolverData>
struct TurbineExternalData : public TurbineBaseData
{
    amrex::Real density{1.0};

    SolverTurbine ext_data;
    ::ext_turb::ExtTurbIface<SolverTurbine, SolverData>* ext_ptr{nullptr};

    MPI_Comm tcomm{MPI_COMM_NULL};
};

template <typename SolverTurbine, typename SolverData>
struct TurbineFromExtSolver : public TurbineType
{
    using InfoType = TurbineInfo;
    using GridType = ActGrid;
    using MetaType = TurbineExternalData<SolverTurbine, SolverData>;
    using DataType =
        ActDataHolder<TurbineFromExtSolver<SolverTurbine, SolverData>>;

    static std::string identifier() { return ext_turb::ext_id<SolverTurbine>(); }
};

} // namespace amr_wind::actuator

#endif /* TURBINEEXTERNAL_H */
