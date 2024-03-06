#ifndef ShearStress_H
#define ShearStress_H

#include "amr-wind/wind_energy/MOData.H"

/**
 *  \defgroup shear_stress Shear Stress Wall Models
 *
 *   ShearStress contains functions to compute velocity and temperature shear
 * stress wall models the default is the Moeng wall model specifying the wall
 * model is done through the input file using ABL.wall_shear_stress_type options
 * include "constant", "local", "Schumann", and "Moeng"
 *
 * \ingroup we_abl
 */

namespace amr_wind {

struct ShearStressConstant
{
    explicit ShearStressConstant(const amr_wind::MOData& mo)
        : utau2(mo.utau * mo.utau)
        , u_mean(mo.vel_mean[0])
        , v_mean(mo.vel_mean[1])
        , wspd_mean(mo.vmag_mean)
        , theta_mean(mo.theta_mean)
        , theta_surface(mo.surf_temp)
        , term1(mo.utau * mo.kappa / mo.phi_h())
    {}

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
    calc_vel_x(amrex::Real /* u */, amrex::Real /* wspd */) const
    {
        return u_mean / wspd_mean * utau2;
    };

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
    calc_vel_y(amrex::Real /* u */, amrex::Real /* wspd */) const
    {
        return v_mean / wspd_mean * utau2;
    };

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
    calc_theta(amrex::Real /* wspd */, amrex::Real /* theta */) const
    {
        return term1 * (theta_mean - theta_surface);
    };

    amrex::Real utau2;
    amrex::Real u_mean;
    amrex::Real v_mean;
    amrex::Real wspd_mean;
    amrex::Real theta_mean;
    amrex::Real theta_surface;
    amrex::Real term1;
};

struct ShearStressLocal
{
    explicit ShearStressLocal(const amr_wind::MOData& mo)
        : utau2(mo.utau * mo.utau)
        , theta_surface(mo.surf_temp)
        , term1(mo.utau * mo.kappa / mo.phi_h())
    {}

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
    calc_vel_x(amrex::Real u, amrex::Real wspd) const
    {
        return u / amrex::max(wspd, small_vel) * utau2;
    };

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
    calc_vel_y(amrex::Real v, amrex::Real wspd) const
    {
        return v / amrex::max(wspd, small_vel) * utau2;
    };

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
    calc_theta(amrex::Real /* wspd */, amrex::Real theta) const
    {
        return term1 * (theta - theta_surface);
    };

    amrex::Real utau2;
    amrex::Real theta_surface;
    amrex::Real term1;
    amrex::Real small_vel{1.0e-6};
};

struct ShearStressSchumann
{
    explicit ShearStressSchumann(const amr_wind::MOData& mo)
        : utau2(mo.utau * mo.utau)
        , wspd_mean(mo.vmag_mean)
        , theta_surface(mo.surf_temp)
        , term1(mo.utau * mo.kappa / mo.phi_h())
    {}

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
    calc_vel_x(amrex::Real u, amrex::Real /* wspd */) const
    {
        return u / wspd_mean * utau2;
    };

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
    calc_vel_y(amrex::Real v, amrex::Real /* wspd */) const
    {
        return v / wspd_mean * utau2;
    };

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
    calc_theta(amrex::Real /* wspd */, amrex::Real theta) const
    {
        return term1 * (theta - theta_surface);
    };

    amrex::Real utau2;
    amrex::Real wspd_mean;
    amrex::Real theta_surface;
    amrex::Real term1;
};

struct ShearStressMoeng
{
    explicit ShearStressMoeng(const amr_wind::MOData& mo)
        : utau2(mo.utau * mo.utau)
        , u_mean(mo.vel_mean[0])
        , v_mean(mo.vel_mean[1])
        , wspd_mean(mo.vmag_mean)
        , theta_surface(mo.surf_temp)
        , theta_mean(mo.theta_mean)
        , term1(mo.utau * mo.kappa / (mo.vmag_mean * mo.phi_h()))
    {}

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
    calc_vel_x(amrex::Real u, amrex::Real wspd) const
    {
        return ((u - u_mean) * wspd_mean + wspd * u_mean) /
               (wspd_mean * wspd_mean) * utau2;
    };

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
    calc_vel_y(amrex::Real v, amrex::Real wspd) const
    {
        return ((v - v_mean) * wspd_mean + wspd * v_mean) /
               (wspd_mean * wspd_mean) * utau2;
    };

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
    calc_theta(amrex::Real wspd, amrex::Real theta) const
    {
        const amrex::Real num1 = (theta - theta_mean) * wspd_mean;
        const amrex::Real num2 = (theta_mean - theta_surface) * wspd;
        return term1 * (num1 + num2);
    };

    amrex::Real utau2;
    amrex::Real u_mean;
    amrex::Real v_mean;
    amrex::Real wspd_mean;
    amrex::Real theta_surface;
    amrex::Real theta_mean;
    amrex::Real term1;
};

struct ShearStressDonelan
{
    explicit ShearStressDonelan(const amr_wind::MOData& mo)
        : wspd_mean(mo.vmag_mean)
        , theta_mean(mo.theta_mean)
        , theta_surface(mo.surf_temp)
        , temp_flux_surface(mo.surf_temp_flux)
        , alg_type(mo.alg_type)
    {}

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
    calc_vel_x(amrex::Real u, amrex::Real wspd) const
    {
        amrex::Real Cd = 0.001;
        const amrex::Real c = 7e-5;
        if (wspd_mean <= 5.0) {
            Cd = 0.001;
        } else if (wspd_mean < 25.0 && wspd_mean > 5.0) {
            Cd = 0.001 + c * (wspd_mean - 5.0);
        } else {
            Cd = 0.0024;
        }

        return Cd * u * wspd;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
    calc_vel_y(amrex::Real v, amrex::Real wspd) const
    {
        amrex::Real Cd = 0.001;
        const amrex::Real c = 7e-5;
        if (wspd_mean <= 5) {
            Cd = 0.001;
        } else if (wspd_mean < 25.0 && wspd_mean > 5.0) {
            Cd = 0.001 + c * (wspd_mean - 5);
        } else {
            Cd = 0.0024;
        }

        return Cd * v * wspd;
    };

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
    calc_theta(amrex::Real /* wspd */, amrex::Real /* theta */) const
    {
        amrex::Real flux = 0.0;
        switch (alg_type) {
        case amr_wind::MOData::ThetaCalcType::HEAT_FLUX:
            flux = temp_flux_surface;
            break;
        case amr_wind::MOData::ThetaCalcType::SURFACE_TEMPERATURE:
            flux = 0.0012 * wspd_mean * (theta_surface - theta_mean);
            break;
        }
        return flux;
    };

    amrex::Real wspd_mean;
    amrex::Real theta_mean;
    amrex::Real theta_surface;
    amrex::Real temp_flux_surface;
    amr_wind::MOData::ThetaCalcType alg_type;
};

} // namespace amr_wind

#endif /* ShearStress_H */
