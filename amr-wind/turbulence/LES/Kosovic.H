#ifndef KOSOVIC_H
#define KOSOVIC_H

#include <string>
#include "amr-wind/turbulence/TurbModelBase.H"
#include "amr-wind/core/FieldRepo.H"

namespace amr_wind::turbulence {

/** Kosovic LES Model
 *  \ingroup turb_model
 */
template <typename Transport>
class Kosovic : public TurbModelBase<Transport>
{
public:
    static std::string identifier()
    {
        return "Kosovic-" + Transport::identifier();
    }

    explicit Kosovic(CFDSim& sim);

    //! Model name for debugging purposes
    std::string model_name() const override { return "Kosovic"; }

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
    //! Kosovic coefficient (default value set for ABL simulations)
  amrex::Real m_Cb{0.36};
  amrex::Real m_Cs{0.135};
  amrex::Real m_C1{2.1};
  amrex::Real m_C2{2.1};
  amrex::Real m_Sk{0.5};
  bool m_surfaceRANS{false};
  amrex::Real m_surfaceFactor{0};
  amrex::Real m_switchLoc{24};
  amrex::Real m_surfaceRANSExp{2};
  amrex::Real m_LESTurnOff{1e15};
  bool m_writeTerms{false};
  amrex::Real m_refMOL{-1e30};
  const Field& m_vel;
  const Field& m_rho;
  Field& m_Nij;
  Field& m_divNij;
};

} // namespace amr_wind::turbulence

#endif /* SMAGORINSKY_H */
