#ifndef KLAXELL_H
#define KLAXELL_H

#include <string>
#include "amr-wind/turbulence/TurbModelBase.H"

namespace amr_wind::turbulence {

/** Base class for 1-Equation RANS TKE turbulence model
 *  \ingroup turb_model
 */
template <typename Transport>
class KLAxell : public TurbModelBase<Transport>
{
public:
    static std::string identifier()
    {
        return "KLAxell-" + Transport::identifier();
    }

    explicit KLAxell(CFDSim& sim);

    std::string model_name() const override { return "KLAxell"; }

    //! Update the turbulent viscosity field
    void update_turbulent_viscosity(
        const FieldState fstate, const DiffusionType /*unused*/) override;

    //! Do any post advance work
    void post_advance_work() override;

    //! Update the effective thermal diffusivity field
    void update_alphaeff(Field& alphaeff) override;

    //! Update the effective scalar diffusivity field
    void update_scalar_diff(Field& deff, const std::string& name) override;

    //! Parse turbulence model coefficients
    void parse_model_coeffs() override;

    //! Return turbulence model coefficients
    TurbulenceModel::CoeffsDictType model_coeffs() const override;

private:
    Field& m_vel;
    Field& m_turb_lscale;
    Field& m_shear_prod;
    Field& m_buoy_prod;
    Field& m_dissip;
    Field& m_rho;
    Field* m_tke{nullptr};
    //! Turbulence constant
    amrex::Real m_Cmu{0.556};
    amrex::Real m_Cmu_prime{0.556};
    amrex::Real m_prandtl{1.0};
    Field& m_temperature;
    //! Gravity vector (m/s^2)
    amrex::Vector<amrex::Real> m_gravity{0.0, 0.0, -9.81};
    //! Reference temperature (Kelvin)
    amrex::Real m_ref_theta{300.0};
    amrex::Real m_surf_flux{0};
    amrex::Real m_lengthscale_switch{800};
};

} // namespace amr_wind::turbulence

#endif /* KLAXELL_H */
