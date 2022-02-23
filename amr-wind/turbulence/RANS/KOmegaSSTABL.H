#ifndef KOMEGASSTABL_H
#define KOMEGASSTABL_H

#include <string>
#include "amr-wind/turbulence/RANS/KOmegaSST.H"

namespace amr_wind {
namespace turbulence {

/** K-Omega-SST RANS model for ABL
 *
 *
 *
 */
template <typename Transport>
class KOmegaSSTABL : public KOmegaSST<Transport>
{
public:
    static std::string identifier()
    {
        return "KOmegaSSTABL-" + Transport::identifier();
    }

    explicit KOmegaSSTABL(CFDSim& sim);

    virtual ~KOmegaSSTABL();

    virtual std::string model_name() const override { return "KOmegaSSTABL"; }

    //! Update the turbulent viscosity field
    virtual void update_turbulent_viscosity(const FieldState fstate) override;

    //! No post advance work for this model
    virtual void post_advance_work() override {}

    //! Update the effective scalar diffusivity field
    virtual void
    update_scalar_diff(Field& deff, const std::string& name) override;

    //! Return turbulence model coefficients
    TurbulenceModel::CoeffsDictType model_coeffs() const override;

protected:
    //! Turbulence constants
    //  Add extra vars in here
};

} // namespace turbulence
} // namespace amr_wind

#endif /* KOMEGASSTABL_H */
