#ifndef TWOPHASETRANSPORT_H
#define TWOPHASETRANSPORT_H

#include "amr-wind/transport_models/TransportModel.H"
#include "AMReX_ParmParse.H"

namespace amr_wind {
namespace transport {

/** Two-phase fluid transport model with piece-wise constant properties
 *  \ingroup transport
 */
class TwoPhaseTransport : public TransportModel
{
public:
    static constexpr bool constant_properties = false;

    static const std::string identifier() { return "TwoPhaseTransport"; }

    TwoPhaseTransport(FieldRepo& repo) : m_repo(repo)
    {
        amrex::ParmParse pp("transport");
        pp.query("viscosity_fluid1", m_mu1);
        pp.query("viscosity_fluid2", m_mu2);
        pp.query("laminar_prandtl", m_Pr);
        pp.query("turbulent_prandtl", m_Prt);
    }

    virtual ~TwoPhaseTransport() = default;

    inline amrex::Real laminar_prandtl() const { return m_Pr; }

    inline amrex::Real turbulent_prandtl() const { return m_Prt; }

    inline amrex::Real laminar_schmidt(const std::string& scalar_name) const
    {
        amrex::ParmParse pp("transport");
        const std::string key = scalar_name + "_laminar_schmidt";
        amrex::Real lam_schmidt = 1.0;
        pp.query(key.c_str(), lam_schmidt);
        return lam_schmidt;
    }

    inline amrex::Real turbulent_schmidt(const std::string& scalar_name) const
    {
        amrex::ParmParse pp("transport");
        const std::string key = scalar_name + "_turbulent_schmidt";
        amrex::Real turb_schmidt = 1.0;
        pp.query(key.c_str(), turb_schmidt);
        return turb_schmidt;
    }

    //! Return the dynamic visocity field
    inline std::unique_ptr<ScratchField> mu()
    {
        auto mu = m_repo.create_scratch_field(1, 1);
        auto& levelset = m_repo.get_field("levelset");
        const auto& geom=m_repo.mesh().Geom();

        for (int lev = 0; lev < m_repo.num_active_levels(); ++lev) {
            for (amrex::MFIter mfi((*mu)(lev)); mfi.isValid(); ++mfi) {
                const auto& vbx = mfi.validbox();
                const auto& dx = geom[lev].CellSizeArray();

                const amrex::Array4<amrex::Real>& phi = levelset(lev).array(mfi);
                const amrex::Array4<amrex::Real>& visc = (*mu)(lev).array(mfi);
                const amrex::Real eps = std::cbrt(2.*dx[0] * dx[1] * dx[2]);
                const amrex::Real mu1 = m_mu1;
                const amrex::Real mu2 = m_mu2;
                amrex::ParallelFor(
                vbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    amrex::Real H;
                    if (phi(i, j, k) > eps) {
                        H = 1.0;
                    } else if (phi(i, j, k) < -eps) {
                        H = 0.;
                    } else {
                        H = 0.5 * (1 + phi(i, j, k) / eps +
                               1. / M_PI * std::sin(phi(i, j, k) * M_PI / eps));
                    }
                    visc(i, j, k) = mu1 * H + mu2 * (1 - H);
            });
            }
        }
        return mu;
    }

    //! Return the thermal diffusivity field
    inline std::unique_ptr<ScratchField> alpha()
    {
        auto alpha = mu();
        amrex::Real inv_Pr = 1.0 / m_Pr;
        for (int lev = 0; lev < m_repo.num_active_levels(); ++lev) {
            (*alpha)(lev).mult(inv_Pr);
        }
        return alpha;
    }

    inline std::unique_ptr<ScratchField>
    scalar_diffusivity(const std::string& scalar_name)
    {
        amrex::Real lam_schmidt = laminar_schmidt(scalar_name);

        amrex::Real inv_schmidt = 1.0 / lam_schmidt;
        auto diff = mu();
        for (int lev = 0; lev < m_repo.num_active_levels(); ++lev) {
            (*diff)(lev).mult(inv_schmidt);
        }

        return diff;
    }

private:
    //! Reference to the field repository (for creating scratch fields)
    FieldRepo& m_repo;

    //! Phase 1 (Laminar) dynamic viscosity
    amrex::Real m_mu1{1.0e-4};

    //! Phase 2 (Laminar) dynamic viscosity
    amrex::Real m_mu2{1.0e-5};
    
    //! Prandtl number
    amrex::Real m_Pr{1.0};

    //! Turbulent Prandtl number
    amrex::Real m_Prt{1.0};
};

} // namespace transport
} // namespace amr_wind

#endif /* TWOPHASETRANSPORT_H */
