#ifndef KOMEGASSTIDDES_H
#define KOMEGASSTIDDES_H

#include <string>
#include "amr-wind/turbulence/RANS/KOmegaSST.H"

namespace amr_wind::turbulence {

/** K-Omega-SST-IDDES Hybrid RANS/LES turbulence model
 *
 *
 * Turbulence model described in
 *
 * Gritskevich, M. S., Garbaruk, A.V., Schutze, J and Menter, F.R., "Development
 * of DDES and IDDES formulations for the k-omega Shear Stress Transport Model",
 * Flow Turbulence and Combustion, 88, 2012, pp. 431-449, DOI
 * 10.1007/s10494-011-9378-4.
 *
 * Implementation contains the simplified version where f_e = 0 as in eq. 12 in
 * the paper
 *
 */
template <typename Transport>
class KOmegaSSTIDDES : public KOmegaSST<Transport>
{
public:
    static std::string identifier()
    {
        return "KOmegaSSTIDDES-" + Transport::identifier();
    }

    explicit KOmegaSSTIDDES(CFDSim& sim);

    ~KOmegaSSTIDDES() override;

    std::string model_name() const override { return "KOmegaSSTIDDES"; }

    //! Update the turbulent viscosity field
    void update_turbulent_viscosity(const FieldState fstate) override;

    //! No post advance work for this model
    void post_advance_work() override {}

    //! Parse turbulence model coefficients
    void parse_model_coeffs() override;

    //! Return turbulence model coefficients
    TurbulenceModel::CoeffsDictType model_coeffs() const override;

protected:
    //! Turbulence constants
    amrex::Real m_Cdes1{0.78};
    amrex::Real m_Cdes2{0.61};
    amrex::Real m_Cdt1{20.0};
    amrex::Real m_Cdt2{3.0};
    amrex::Real m_Cw{0.15};
    amrex::Real m_kappa{0.41};
};

} // namespace amr_wind::turbulence

#endif /* KOMEGASSTIDDES_H */
