#ifndef ONEEQKSGS_H
#define ONEEQKSGS_H

#include <string>
#include "amr-wind/turbulence/TurbModelBase.H"

namespace amr_wind {
namespace turbulence {

/** Base class for 1-Equation subgrid scale TKE turbulence model
 */
template <typename Transport>
class OneEqKsgs : public TurbModelBase<Transport>
{
public:
    OneEqKsgs(CFDSim& sim);

    virtual ~OneEqKsgs();

protected:
    Field& m_vel;

    Field& m_turb_lscale;
    Field& m_shear_prod;
    Field& m_buoy_prod;
    Field& m_rho;

    Field* m_tke;

    //! Turbulence constant
    amrex::Real m_Ce{0.1};
    amrex::Real m_Ceps{0.93};

};

template <typename Transport>
class OneEqKsgsM84 : public OneEqKsgs<Transport>
{
public:
    static std::string identifier()
    {
        return "OneEqKsgsM84-" + Transport::identifier();
    }

    OneEqKsgsM84(CFDSim& sim);

    virtual ~OneEqKsgsM84();

    virtual std::string model_name() const override { return "OneEqKsgsM84"; }

    //! Update the turbulent viscosity field
    virtual void update_turbulent_viscosity(const FieldState fstate) override;

    //! Update the effective thermal diffusivity field
    virtual void update_alphaeff(Field& alphaeff) override;
    
    //! Return turbulence model coefficients
    TurbulenceModel::CoeffsDictType model_coeffs() const override;

private:

    Field& m_temperature;
    
    //! Gravity vector (m/s^2)
    amrex::Vector<amrex::Real> m_gravity{{0.0, 0.0, -9.81}};

    //! Reference temperature (Kelvin)
    amrex::Real m_ref_theta{300.0};
    
};

template <typename Transport>
class OneEqKsgsS94 : public OneEqKsgs<Transport>
{
public:
    static std::string identifier()
    {
        return "OneEqKsgsS94-" + Transport::identifier();
    }

    OneEqKsgsS94(CFDSim& sim);

    virtual ~OneEqKsgsS94();

    virtual std::string model_name() const override { return "OneEqKsgsS94"; }

    //! Update the turbulent viscosity field
    virtual void update_turbulent_viscosity(const FieldState fstate) override;

    //! Return turbulence model coefficients
    TurbulenceModel::CoeffsDictType model_coeffs() const override;
};

} // namespace turbulence
} // namespace amr_wind

#endif /* ONEEQKSGS_H */
