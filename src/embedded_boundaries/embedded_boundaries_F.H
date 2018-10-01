#include <AMReX_REAL.H>

#include <AMReX_BLFort.H>
#include <AMReX_SPACE.H>

#ifdef BL_SINGLE_PRECISION_PARTICLES
typedef float amrex_particle_real;
#else
typedef double amrex_particle_real;
#endif

#ifdef __cplusplus
extern "C"
{
#endif

    void incflo_get_walls(
            int* bcv, int* exists,
            void* normal, void* center
        );

    void incflo_get_real_walls(
            int* bcv, int* exists,
            void* normal, void* center
        );

    void incflo_eb_to_polygon(
            const amrex::Real* dx,
            const int* slo, const int* shi,
            const void* flag,        const int* fglo, const int* fghi,
            const amrex_real* bcent, const int* blo,  const int* bhi,
            const amrex_real* ax,    const int* axlo, const int* axhi,
            const amrex_real* ay,    const int* aylo, const int* ayhi,
            const amrex_real* az,    const int* azlo, const int* azhi
        );

    void incflo_write_eb_vtp(int* myID);
    void incflo_write_pvtp(int* nProcs);

    void incflo_eb_grid_coverage(
            int* myID, const amrex::Real* dx,
            const int* slo, const int* shi,
            const void* flag, const int* fglo, const int* fghi
        );

#ifdef __cplusplus
}
#endif
