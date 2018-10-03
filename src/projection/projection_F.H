#include <AMReX_REAL.H>
#include <AMReX_BLFort.H>
#include <AMReX_SPACE.H>

#ifdef __cplusplus
extern "C"
{
#endif
   void set_ppe_bc(int* bc_lo, int* bc_hi,
                   const int* domlo, const int* domhi,
                   const int* ng, 
                   const int* bct_ilo, const int* bct_ihi, 
                   const int* bct_jlo, const int* bct_jhi,
                   const int* bct_klo, const int* bct_khi);
#ifdef __cplusplus
}
#endif
