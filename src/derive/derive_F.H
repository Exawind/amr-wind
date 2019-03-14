#include <AMReX_REAL.H>
#include <AMReX_BLFort.H>
#include <AMReX_SPACE.H>

#ifdef __cplusplus
extern "C"
{
#endif
    void compute_strainrate(const int* lo, const int* hi,
                            amrex::Real* sr, const int* slo, const int* shi,
                            amrex::Real* vel, const int* ulo, const int* uhi,
                            const amrex::Real* dx);

    void compute_strainrate_eb(const int* lo, const int* hi,
                               amrex::Real* sr, const int* slo, const int* shi,
                               amrex::Real* vel, const int* ulo, const int* uhi,
                               const void* flag, const int* fglo, const int* fghi,
                               const amrex::Real* dx);

    void compute_vort(const int* lo, const int* hi,
                      amrex::Real* vort, const int* slo, const int* shi,
                      amrex::Real* vel, const int* ulo, const int* uhi,
                      const amrex::Real* dx);

    void compute_vort_eb(const int* lo, const int* hi,
                         amrex::Real* vort, const int* slo, const int* shi,
                         amrex::Real* vel, const int* ulo, const int* uhi,
                         const void* flag, const int* fglo, const int* fghi,
                         const amrex::Real* dx);

    void state_error(const int* lo, const int* hi,
                     char* tag, const int* tag_lo, const int* tag_hi,
                     const BL_FORT_FAB_ARG_3D(state),
                     const char* tagval, const char* clearval,
                     const amrex_real* dx, const amrex_real* problo, const amrex_real* time);
#ifdef __cplusplus
}
#endif
