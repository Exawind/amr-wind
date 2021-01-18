#ifndef KOMEGASST_H
#define KOMEGASST_H

#include <string>
#include "amr-wind/turbulence/TurbModelBase.H"

namespace amr_wind {
namespace turbulence {

/** Base class for 2-Equation k-w RANS turbulence model
 */
template <typename Transport>
class KOmegaSST : public TurbModelBase<Transport>
{
public:
    static std::string identifier()
    {
        return "KOmegaSST-" + Transport::identifier();
    }
    
    KOmegaSST(CFDSim& sim);

    virtual ~KOmegaSST();

    virtual std::string model_name() const override { return "KOmegaSST"; }

    //! Update the turbulent viscosity field
    virtual void update_turbulent_viscosity(const FieldState fstate) override;

    //! Update the effective thermal diffusivity field
    virtual void update_alphaeff(Field& alphaeff) override;
    
    //! Return turbulence model coefficients
    TurbulenceModel::CoeffsDictType model_coeffs() const override;
    
protected:
    Field& m_vel;

    Field& m_shear_prod;
    Field& m_diss;
    Field& m_rho;

    Field* m_tke;
    Field* m_omega;
    Field* m_walldist;

    //! Turbulence constants
    amrex::Real m_beta_star{0.09};
    amrex::Real m_alpha1{0.5555555555555556};
    amrex::Real m_alpha2{0.44};
    amrex::Real m_beta1{0.075};
    amrex::Real m_beta2{0.0828};
    amrex::Real m_sigma_k1{0.85};
    amrex::Real m_sigma_k2{1.0};
    amrex::Real m_sigma_omega1{0.5};
    amrex::Real m_sigma_omega2{0.856};

};


} // namespace turbulence
} // namespace amr_wind

#endif /* KOMEGASST_H */
