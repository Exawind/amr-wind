#ifndef SMAGORINSKY_H
#define SMAGORINSKY_H

#include <string>
#include "amr-wind/turbulence/TurbModelBase.H"
#include "amr-wind/core/FieldRepo.H"

namespace amr_wind::turbulence {

/** Smagorinsky LES Model
 *  \ingroup turb_model
 */
template <typename Transport>
class Smagorinsky : public TurbModelBase<Transport>
{
public:
    static std::string identifier()
    {
        return "Smagorinsky-" + Transport::identifier();
    }

    explicit Smagorinsky(CFDSim& sim);

    //! Model name for debugging purposes
    std::string model_name() const override { return "Smagorinsky"; }

    //! Update the turbulent viscosity field
    void update_turbulent_viscosity(
        const FieldState fstate, const DiffusionType /*unused*/) override;

    //! No post advance work for this model
    void post_advance_work() override {}

    //! Parse turbulence model coefficients
    void parse_model_coeffs() override;

    //! Return model coefficients dictionary
    TurbulenceModel::CoeffsDictType model_coeffs() const override;

private:
    //! Smagorinsky coefficient (default value set for ABL simulations)
    amrex::Real m_Cs{0.135};

    const Field& m_vel;
    const Field& m_rho;
};

} // namespace amr_wind::turbulence

#endif /* SMAGORINSKY_H */
