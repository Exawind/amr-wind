//
//  ShearStress.H
//  amr-wind
//
//

#ifndef ShearStress_H
#define ShearStress_H

#include "amr-wind/wind_energy/MOData.H"

struct ShearStressConstant
{
    ShearStressConstant(amr_wind::MOData mo)
        : utau2(mo.utau * mo.utau)
        , u_mean(mo.vel_mean[0])
        , v_mean(mo.vel_mean[1])
        , wspd_mean(mo.vmag_mean)
        , theta_mean(mo.theta_mean)
        , theta_surface(mo.surf_temp)
    {
        term1 = mo.utau * mo.kappa / mo.phi_h();
    }

    AMREX_GPU_HOST_DEVICE amrex::Real
        calc_vel_x(amrex::Real /* u */, amrex::Real /* wspd */) const
    {
        return u_mean / wspd_mean * utau2;
    };

    AMREX_GPU_HOST_DEVICE amrex::Real
        calc_vel_y(amrex::Real /* u */, amrex::Real /* wspd */) const
    {
        return v_mean / wspd_mean * utau2;
    };

    AMREX_GPU_HOST_DEVICE amrex::Real
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
    ShearStressLocal(amr_wind::MOData mo)
        : utau2(mo.utau * mo.utau), theta_surface(mo.surf_temp)
    {
        term1 = mo.utau * mo.kappa / mo.phi_h();
    }

    AMREX_GPU_HOST_DEVICE amrex::Real
    calc_vel_x(amrex::Real u, amrex::Real wspd) const
    {
        return u / amrex::max(wspd, 1.0e-6) * utau2;
    };

    AMREX_GPU_HOST_DEVICE amrex::Real
    calc_vel_y(amrex::Real v, amrex::Real wspd) const
    {
        return v / amrex::max(wspd, 1.0e-6) * utau2;
    };

    AMREX_GPU_HOST_DEVICE amrex::Real
    calc_theta(amrex::Real /* wspd */, amrex::Real theta) const
    {
        return term1 * (theta - theta_surface);
    };

    amrex::Real utau2;
    amrex::Real theta_surface;
    amrex::Real term1;
};

struct ShearStressSchumann
{
    ShearStressSchumann(amr_wind::MOData mo)
        : utau2(mo.utau * mo.utau)
        , wspd_mean(mo.vmag_mean)
        , theta_surface(mo.surf_temp)
    {
        term1 = mo.utau * mo.kappa / mo.phi_h();
    }

    AMREX_GPU_HOST_DEVICE amrex::Real
    calc_vel_x(amrex::Real u, amrex::Real /* wspd */) const
    {
        return u / wspd_mean * utau2;
    };

    AMREX_GPU_HOST_DEVICE amrex::Real
    calc_vel_y(amrex::Real v, amrex::Real /* wspd */) const
    {
        return v / wspd_mean * utau2;
    };

    AMREX_GPU_HOST_DEVICE amrex::Real
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
    ShearStressMoeng(amr_wind::MOData mo)
        : utau2(mo.utau * mo.utau)
        , u_mean(mo.vel_mean[0])
        , v_mean(mo.vel_mean[1])
        , wspd_mean(mo.vmag_mean)
        , theta_surface(mo.surf_temp)
        , theta_mean(mo.theta_mean)
    {
        term1 = mo.utau * mo.kappa / (mo.vmag_mean * mo.phi_h());
    }

    AMREX_GPU_HOST_DEVICE amrex::Real
    calc_vel_x(amrex::Real u, amrex::Real wspd) const
    {
        return ((u - u_mean) * wspd_mean + wspd * u_mean) /
               (wspd_mean * wspd_mean) * utau2;
    };

    AMREX_GPU_HOST_DEVICE amrex::Real
    calc_vel_y(amrex::Real v, amrex::Real wspd) const
    {
        return ((v - v_mean) * wspd_mean + wspd * v_mean) /
               (wspd_mean * wspd_mean) * utau2;
    };

    AMREX_GPU_HOST_DEVICE amrex::Real
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

#endif /* ShearStress_H */
