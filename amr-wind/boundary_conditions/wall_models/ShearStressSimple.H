#ifndef SHEARSTRESSSIMPLE_H
#define SHEARSTRESSSIMPLE_H

#include "amr-wind/boundary_conditions/wall_models/LogLaw.H"
#include "amr-wind/wind_energy/ShearStress.H"

namespace amr_wind {

struct SimpleShearSchumann
{
    explicit SimpleShearSchumann(const amr_wind::LogLaw& ll)
        : utau2(ll.utau_mean * ll.utau_mean), wspd_mean(ll.wspd_mean), m_ll(ll)
    {}

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
    get_shear(amrex::Real u, amrex::Real /* wspd */) const
    {
        return u / wspd_mean * utau2;
    };

    amrex::Real utau2;
    amrex::Real wspd_mean;
    const amr_wind::LogLaw m_ll;
};
struct SimpleShearLogLaw
{
    explicit SimpleShearLogLaw(const amr_wind::LogLaw& ll) : m_ll(ll) {}

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
    get_shear(amrex::Real u, amrex::Real wspd) const
    {
        amrex::Real utau = m_ll.get_utau(wspd);
        return utau * utau * u / wspd;
    };

    const amr_wind::LogLaw m_ll;
};
} // namespace amr_wind

#endif /* SHEARSTRESSSIMPLE_H */
