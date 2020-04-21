#ifndef TURBULENCEMODEL_H
#define TURBULENCEMODEL_H

#include "Factory.H"
#include "FieldDescTypes.H"

namespace amr_wind {

class CFDSim;
class ScratchField;
class Field;

namespace turbulence {

class TurbulenceModel : public Factory<TurbulenceModel, CFDSim&>
{
public:
    static std::string base_identifier() { return "TurbulenceModel"; }

    virtual ~TurbulenceModel() = default;

    virtual void update_turbulent_viscosity(const FieldState fstate) = 0;

    virtual void register_mueff_field(Field& mueff) = 0;

    virtual void register_alphaeff_field(Field& alphaeff) = 0;

    virtual std::string model_name() = 0;

    //! Return the dynamic viscosity (laminar) field
    virtual std::unique_ptr<ScratchField> mu() = 0;

    //! Return the thermal diffusivity (laminar) field for enthalpy/temperature
    virtual std::unique_ptr<ScratchField> alpha() = 0;

    //! Return the scalar diffusivity field
    virtual std::unique_ptr<ScratchField> scalar_diffusivity(const std::string& name) = 0;

    //! Return the turbulent dynamic viscosity field
    virtual Field& mu_turb() = 0;

    //! Return the turbulent dynamic viscosity field
    virtual Field& alpha_turb() = 0;

    //! Return the effective dynamic viscosity field
    virtual Field& mueff() = 0;

    //! Return the effective thermal diffusivity field
    virtual Field& alphaeff() = 0;

    //! Flag indicating whether the model is turbulent
    virtual bool is_turbulent() const { return true; }

    //! Interface to update effective viscosity (mu_eff = mu + mu_t)
    virtual void update_mueff(Field& mueff) = 0;

    //! Interface to update effective thermal diffusivity
    virtual void update_alphaeff(Field& alphaeff) = 0;

    //! Interface to update scalar diffusivity based on Schmidt number
    virtual void update_scalar_diff(Field& deff, const std::string& name) = 0;
};

}
}

#endif /* TURBULENCEMODEL_H */
