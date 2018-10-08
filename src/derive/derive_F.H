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
                               const amrex::Real* afrac_x,const int*  axlo,const int*  axhi, 
                               const amrex::Real* afrac_y,const int*  aylo,const int*  ayhi, 
                               const amrex::Real* afrac_z,const int*  azlo,const int*  azhi, 
                               const amrex::Real* cent_x,const int*   cxlo,const int*  cxhi, 
                               const amrex::Real* cent_y,const int*   cylo,const int*  cyhi, 
                               const amrex::Real* cent_z,const int*   czlo,const int*  czhi,
                               const amrex::Real* vfrac, const int* vflo, const int*  vfhi, 
                               const amrex::Real* bcent, const int* blo, const int*  bhi, 
                               const amrex::Real* dx);

    void compute_vort(const int* lo, const int* hi,
                      amrex::Real* vort, const int* slo, const int* shi,
                      amrex::Real* vel, const int* ulo, const int* uhi,
                      const amrex::Real* dx);

    void compute_vort_eb(const int* lo, const int* hi,
                         amrex::Real* vort, const int* slo, const int* shi,
                         amrex::Real* vel, const int* ulo, const int* uhi,
                         const void* flag, const int* fglo, const int* fghi,
                         const amrex::Real* afrac_x,const int*  axlo,const int*  axhi, 
                         const amrex::Real* afrac_y,const int*  aylo,const int*  ayhi, 
                         const amrex::Real* afrac_z,const int*  azlo,const int*  azhi, 
                         const amrex::Real* cent_x,const int*   cxlo,const int*  cxhi, 
                         const amrex::Real* cent_y,const int*   cylo,const int*  cyhi, 
                         const amrex::Real* cent_z,const int*   czlo,const int*  czhi,
                         const amrex::Real* vfrac, const int* vflo, const int*  vfhi, 
                         const amrex::Real* bcent, const int* blo, const int*  bhi, 
                         const amrex::Real* dx);
#ifdef __cplusplus
}
#endif
