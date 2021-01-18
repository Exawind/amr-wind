#ifndef KOMEGASST_H
#define KOMEGASST_H

#include <string>
#include "amr-wind/turbulence/TurbModelBase.H"

namespace amr_wind {
namespace turbulence {

/** K-Omega-SST RANS turbulence model
 * 
 * This also serves as the base class for all k-omega type RANS models
 * 
 * Turbulence model described in
 * 
 * Menter, F., Esch, T.,
 * "Elements of Industrial Heat Transfer Prediction",
 * 16th Brazilian Congress of Mechanical Engineering (COBEM),
 * Nov. 2001.
 * 
 * with updated coefficients from
 * 
 * Menter, F. R., Kuntz, M., and Langtry, R.,
 * "Ten Years of Industrial Experience with the SST Turbulence Model",
 * Turbulence, Heat and Mass Transfer 4, 2003,
 * pp. 625 - 632.
 * 
 * Also see
 * 
 * http://turbmodels.larc.nasa.gov/sst.html
 * 
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

    //! Update the effective scalar diffusivity field
    virtual void update_scalar_diff(Field& deff, const std::string& name) override;
    
    //! Return turbulence model coefficients
    TurbulenceModel::CoeffsDictType model_coeffs() const override;
    
protected:
    Field& m_vel;

    Field& m_f1;
    Field& m_shear_prod;
    Field& m_diss;    
    Field& m_sdr_src;
    Field& m_sdr_diss;    
    Field& m_rho;
    Field& m_walldist;
    
    Field* m_tke;
    Field* m_sdr;


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
    amrex::Real m_a1{0.31};
 
};

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
class KOmegaSSTIDDES: public KOmegaSST<Transport>
{
public:
    static std::string identifier()
    {
        return "KOmegaSSTIDDES-" + Transport::identifier();
    }

    KOmegaSSTIDDES(CFDSim& sim);

    virtual ~KOmegaSSTIDDES();

    virtual std::string model_name() const override { return "KOmegaSSTIDDES"; }

    //! Update the turbulent viscosity field
    virtual void update_turbulent_viscosity(const FieldState fstate) override;
    
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


} // namespace turbulence
} // namespace amr_wind

#endif /* KOMEGASST_H */
