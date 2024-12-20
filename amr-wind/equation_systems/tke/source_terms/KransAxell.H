#ifndef KRANSAXELL_H
#define KRANSAXELL_H

#include "amr-wind/equation_systems/tke/TKESource.H"

namespace amr_wind::pde::tke {

/** TKE source term based on Axell 2011 paper
 * Axell, L. B., & Liungman, O. (2001). A one-equation turbulence model for
 * geophysical applications: comparison with data and the k− ε model.
 * Environmental Fluid Mechanics, 1, 71-106.
 *  \ingroup tke_src turb_model we_abl
 */
class KransAxell : public TKESource::Register<KransAxell>
{
public:
    static std::string identifier() { return "KransAxell"; }

    explicit KransAxell(const CFDSim& /*sim*/);

    ~KransAxell() override;

    void operator()(
        const int lev,
        const amrex::MFIter& mfi,
        const amrex::Box& bx,
        const FieldState fstate,
        const amrex::Array4<amrex::Real>& src_term) const override;

private:
    Field& m_turb_lscale;
    Field& m_shear_prod;
    Field& m_buoy_prod;
    Field& m_dissip;
    Field& m_tke;
    amrex::Real m_Cmu{0.556};
    amrex::Real m_heat_flux{0.0};
    amrex::Real m_ref_temp{300.0};
    amrex::Real m_z0{0.1};
    amrex::Real m_kappa{0.41};
    amrex::Real m_sponge_start{600};
    amrex::Real m_ref_tke{1e-10};
    amrex::Vector<amrex::Real> m_gravity{0.0, 0.0, -9.81};
    const SimTime& m_time;
    const CFDSim& m_sim;
    const amrex::AmrCore& m_mesh;
    const Field& m_velocity;
    std::string m_wall_het_model{"none"};
    amrex::Real m_mol_length{10000};
    amrex::Real m_gamma_m{5.0};
    amrex::Real m_beta_m{16.0};

    amrex::Real stability(
        const amrex::Real z,
        const amrex::Real mol_length,
        const amrex::Real gamma_m,
        const amrex::Real beta_m) const
    {
        const amrex::Real zeta = z / mol_length;
        amrex::Real psi_m = 0;
        if (zeta > 0) {
            psi_m = -gamma_m * zeta;
        } else {
            amrex::Real x = std::sqrt(std::sqrt(1 - beta_m * zeta));
            psi_m = 2.0 * std::log(0.5 * (1.0 + x)) + log(0.5 * (1 + x * x)) -
                    2.0 * std::atan(x) + 2.0 * std::atan(1.0);
        }
        return psi_m;
    }
};

} // namespace amr_wind::pde::tke

#endif /* KRANSAXELL_H */
