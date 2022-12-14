#ifndef TWOPHASETRANSPORT_H
#define TWOPHASETRANSPORT_H

#include "amr-wind/physics/multiphase/MultiPhase.H"
#include "amr-wind/transport_models/TransportModel.H"
#include "AMReX_ParmParse.H"

namespace amr_wind::transport {

/** Two-phase fluid transport model with piece-wise constant properties
 *  \ingroup transport
 */
class TwoPhaseTransport : public TransportModel
{
public:
    static constexpr bool constant_properties = false;

    static std::string identifier() { return "TwoPhaseTransport"; }

    explicit TwoPhaseTransport(CFDSim& sim) : m_sim(sim), m_repo(sim.repo())
    {
        auto& physics_mgr = m_sim.physics_manager();
        if (!physics_mgr.contains("MultiPhase")) {
            amrex::Abort("TwoPhaseTransport requires MultiPhase physics");
        }
        auto& multiphase = physics_mgr.get<MultiPhase>();
        m_ifacetype = multiphase.interface_capturing_method();

        amrex::ParmParse pp("transport");
        pp.query("viscosity_fluid1", m_mu1);
        pp.query("viscosity_fluid2", m_mu2);
        pp.query("laminar_prandtl", m_Pr);
        pp.query("turbulent_prandtl", m_Prt);
    }

    ~TwoPhaseTransport() override = default;

    inline amrex::Real laminar_prandtl() const { return m_Pr; }

    inline amrex::Real turbulent_prandtl() const { return m_Prt; }

    static inline amrex::Real laminar_schmidt(const std::string& scalar_name)
    {
        amrex::ParmParse pp("transport");
        const std::string key = scalar_name + "_laminar_schmidt";
        amrex::Real lam_schmidt = 1.0;
        pp.query(key.c_str(), lam_schmidt);
        return lam_schmidt;
    }

    static inline amrex::Real turbulent_schmidt(const std::string& scalar_name)
    {
        amrex::ParmParse pp("transport");
        const std::string key = scalar_name + "_turbulent_schmidt";
        amrex::Real turb_schmidt = 1.0;
        pp.query(key.c_str(), turb_schmidt);
        return turb_schmidt;
    }

    //! Return the dynamic visocity field
    inline std::unique_ptr<ScratchField> mu() override
    {
        // Select the interface capturing method
        auto mu = m_repo.create_scratch_field(1, 1);

        if (m_ifacetype == InterfaceCapturingMethod::VOF) {

            // cppcheck-suppress constVariable
            auto& vof = m_repo.get_field("vof");

            for (int lev = 0; lev < m_repo.num_active_levels(); ++lev) {
                for (amrex::MFIter mfi((*mu)(lev)); mfi.isValid(); ++mfi) {
                    const auto& vbx = mfi.growntilebox();
                    const amrex::Array4<amrex::Real>& volfrac =
                        vof(lev).array(mfi);
                    const amrex::Array4<amrex::Real>& visc =
                        (*mu)(lev).array(mfi);
                    const amrex::Real mu1 = m_mu1;
                    const amrex::Real mu2 = m_mu2;
                    amrex::ParallelFor(
                        vbx,
                        [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                            visc(i, j, k) = mu1 * volfrac(i, j, k) +
                                            mu2 * (1.0 - volfrac(i, j, k));
                        });
                }
            }

        } else if (m_ifacetype == InterfaceCapturingMethod::LS) {

            // cppcheck-suppress constVariable
            auto& levelset = m_repo.get_field("levelset");
            const auto& geom = m_repo.mesh().Geom();

            for (int lev = 0; lev < m_repo.num_active_levels(); ++lev) {
                for (amrex::MFIter mfi((*mu)(lev)); mfi.isValid(); ++mfi) {
                    const auto& vbx = mfi.growntilebox();
                    const auto& dx = geom[lev].CellSizeArray();
                    const amrex::Array4<amrex::Real>& visc =
                        (*mu)(lev).array(mfi);
                    const amrex::Array4<amrex::Real>& phi =
                        levelset(lev).array(mfi);
                    const amrex::Real eps =
                        std::cbrt(2. * dx[0] * dx[1] * dx[2]);
                    const amrex::Real mu1 = m_mu1;
                    const amrex::Real mu2 = m_mu2;
                    amrex::ParallelFor(
                        vbx,
                        [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                            amrex::Real smooth_heaviside;
                            if (phi(i, j, k) > eps) {
                                smooth_heaviside = 1.0;
                            } else if (phi(i, j, k) < -eps) {
                                smooth_heaviside = 0.;
                            } else {
                                smooth_heaviside =
                                    0.5 *
                                    (1.0 + phi(i, j, k) / eps +
                                     1.0 / M_PI *
                                         std::sin(phi(i, j, k) * M_PI / eps));
                            }
                            visc(i, j, k) = mu1 * smooth_heaviside +
                                            mu2 * (1.0 - smooth_heaviside);
                        });
                }
            }
        }
        return mu;
    }

    //! Return the thermal diffusivity field
    inline std::unique_ptr<ScratchField> alpha() override
    {
        auto alpha = mu();
        amrex::Real inv_Pr = 1.0 / m_Pr;
        for (int lev = 0; lev < m_repo.num_active_levels(); ++lev) {
            (*alpha)(lev).mult(inv_Pr);
        }
        return alpha;
    }

    inline std::unique_ptr<ScratchField>
    scalar_diffusivity(const std::string& scalar_name) override
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
    //! Reference to the CFD sim
    CFDSim& m_sim;

    //! Reference to the field repository (for creating scratch fields)
    FieldRepo& m_repo;

    //! Interface capturing method variable
    InterfaceCapturingMethod m_ifacetype;

    //! Phase 1 (Laminar) dynamic viscosity
    amrex::Real m_mu1{1.0e-4};

    //! Phase 2 (Laminar) dynamic viscosity
    amrex::Real m_mu2{1.0e-5};

    //! Prandtl number
    amrex::Real m_Pr{1.0};

    //! Turbulent Prandtl number
    amrex::Real m_Prt{1.0};
};

} // namespace amr_wind::transport

#endif /* TWOPHASETRANSPORT_H */
