#include <AMReX_REAL.H>

#include <AMReX_BLFort.H>
#include <AMReX_SPACE.H>
#include "AMReX_RealVect.H"

#ifdef __cplusplus
extern "C"
{
#endif
    void check_inputs(amrex::Real* dt);

    void check_initial_conditions(amrex::Real* dx, amrex::Real* dy, amrex::Real* dz,
          const int* domlo, const int* domhi);

    void check_boundary_conditions(amrex::Real* dx, amrex::Real* dy, amrex::Real* dz,
           amrex::Real* xlength, amrex::Real* ylength, amrex::Real* zlength,
           const int* domlo, const int* domhi);

    void check_bc_flow();
#ifdef __cplusplus
}
#endif
