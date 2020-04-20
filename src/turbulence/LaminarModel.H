#ifndef LAMINARMODEL_H
#define LAMINARMODEL_H

#include "TurbModel.H"

namespace amr_wind {
namespace turbulence {

template<typename Transport>
class Laminar : public TurbModel<Transport>
{
public:
    static std::string identifier() { return "Laminar-" + Transport::identifier(); }

    Laminar(CFDSim& sim)
        : TurbModel<Transport>(sim)
    {}

    virtual ~Laminar() = default;

    virtual std::string model_name() { return "laminar"; }

    //! Update the effective/turbulent viscosity field
    virtual void update_turbulent_viscosity() override;

    //! Return the turbulent viscosity field (just the same as the effective field)
    virtual Field& mu_turb() override { return this->mueff(); }

    //! Return the thermal diffusivity field (just the same as effective field)
    virtual Field& alpha_turb() override { return this->alphaeff(); }

    //! Indicate that this model is not a turbulent model type
    virtual bool is_turbulent() const override { return false; }
};

}
}

#endif /* LAMINARMODEL_H */
