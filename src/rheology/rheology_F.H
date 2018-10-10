#include <AMReX_REAL.H>
#include <AMReX_BLFort.H>
#include <AMReX_SPACE.H>

#ifdef __cplusplus
extern "C"
{
#endif
    void powerlaw_viscosity(const int* lo, const int* hi,
                           amrex::Real* eta, const int* elo, const int* ehi,
                           amrex::Real* sr, const int* slo, const int* shi);

    void bingham_viscosity(const int* lo, const int* hi,
                           amrex::Real* eta, const int* elo, const int* ehi,
                           amrex::Real* sr, const int* slo, const int* shi);

    void hb_viscosity(const int* lo, const int* hi,
                           amrex::Real* eta, const int* elo, const int* ehi,
                           amrex::Real* sr, const int* slo, const int* shi);

    void smd_viscosity(const int* lo, const int* hi,
                           amrex::Real* eta, const int* elo, const int* ehi,
                           amrex::Real* sr, const int* slo, const int* shi);
#ifdef __cplusplus
}
#endif
