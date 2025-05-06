#ifndef MULTIPHASEKOSOVIC_H
#define MULTIPHASEKOSOVIC_H

#include <string>
#include "amr-wind/turbulence/TurbModelBase.H"
#include "amr-wind/core/FieldRepo.H"

namespace amr_wind::turbulence {

/** MultiPhaseKosovic LES Model
 *  \ingroup turb_model
 */
template <typename Transport>
class MultiPhaseKosovic : public TurbModelBase<Transport>
{
public:
    static std::string identifier()
    {
        return "MultiPhaseKosovic-" + Transport::identifier();
    }

    explicit MultiPhaseKosovic(CFDSim& sim);

    //! Model name for debugging purposes
    std::string model_name() const override { return "MultiPhaseKosovic"; }

    //! Update the turbulent viscosity field
    void update_turbulent_viscosity(
        const FieldState fstate, const DiffusionType /*unused*/) override;

    //! No post advance work for this model
    void post_advance_work() override {}

    //! Update the effective thermal diffusivity field
    void update_alphaeff(Field& alphaeff) override;

    //! Parse turbulence model coefficients
    void parse_model_coeffs() override;

    //! Return model coefficients dictionary
    TurbulenceModel::CoeffsDictType model_coeffs() const override;

private:
    //! MultiPhaseKosovic coefficient (default value set for ABL simulations)
    // Ref: Mirocha et. al "Implementation of a Nonlinear Subfilter Turbulence
    //  Stress Model for Large-Eddy Simulation in the Advanced Research WRF
    //  Model"
    // , MWR 2012.
    amrex::Real m_Cb{0.36};
    amrex::Real m_Cs{0.135};
    amrex::Real m_C1{2.1};
    amrex::Real m_C2{2.1};
    amrex::Real m_Sk{0.5};
    bool m_writeTerms{false};
    amrex::Vector<amrex::Real> m_gravity{0.0, 0.0, -9.81};
    const Field& m_vel;
    const Field& m_rho;
    Field& m_Nij;
    Field& m_divNij;
};

} // namespace amr_wind::turbulence

#endif /* MULTIPHASEKOSOVIC_H */
