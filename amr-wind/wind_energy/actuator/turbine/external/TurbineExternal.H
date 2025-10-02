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
    ::exw_fast::ExtTurbIface<SolverTurbine,SolverData>* ext_ptr{nullptr};

    MPI_Comm tcomm{MPI_COMM_NULL};
};

template <typename SolverTurbine, typename SolverData>
struct TurbineFromExtSolver : public TurbineType
{
    using InfoType = TurbineInfo;
    using GridType = ActGrid;
    using MetaType = TurbineExternalData<SolverTurbine,SolverData>;
    using DataType = ActDataHolder<TurbineFromExtSolver<SolverTurbine,SolverData>>;

    // !!!! update how identifier is set !!!!!
    static std::string identifier() { return "TurbineFromExtSolver"; }
};

} // namespace amr_wind::actuator

#endif /* TURBINEEXTERNAL_H */
