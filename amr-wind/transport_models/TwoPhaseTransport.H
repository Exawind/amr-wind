#ifndef TWOPHASETRANSPORT_H
#define TWOPHASETRANSPORT_H

#include "amr-wind/physics/multiphase/MultiPhase.H"
#include "amr-wind/transport_models/TransportModel.H"
#include "amr-wind/utilities/constants.H"
#include "AMReX_ParmParse.H"

namespace amr_wind::transport {

/** Two-phase fluid transport model with piece-wise constant properties
 *  \ingroup transport
 */
class TwoPhaseTransport : public TransportModel::Register<TwoPhaseTransport>
{
public:
    static constexpr bool constant_properties = false;

    static std::string identifier() { return "TwoPhaseTransport"; }

    explicit TwoPhaseTransport(const CFDSim& sim)
        : m_sim(sim), m_repo(sim.repo())
    {
        const auto& physics_mgr = m_sim.physics_manager();
        if (!physics_mgr.contains("MultiPhase")) {
            amrex::Abort("TwoPhaseTransport requires MultiPhase physics");
        }
        const auto& multiphase = physics_mgr.get<MultiPhase>();
        m_ifacetype = multiphase.interface_capturing_method();

        amrex::ParmParse pp("transport");
        pp.query("viscosity_fluid1", m_mu1);
        pp.query("viscosity_fluid2", m_mu2);
        pp.query("laminar_prandtl_fluid1", m_Pr1);
        pp.query("laminar_prandtl_fluid2", m_Pr2);
        pp.query("turbulent_prandtl", m_Prt);

        // Check for single-phase quantities and warn
        if (pp.contains("viscosity")) {
            amrex::Print()
                << "WARNING: single-phase viscosity has been specified but "
                   "will not be used! (TwoPhaseTransport)\n";
        }
        if (pp.contains("laminar_prandtl")) {
            amrex::Print()
                << "WARNING: single-phase laminar_prandtl has been specified "
                   "but will not be used! (TwoPhaseTransport)\n";
        }

        // Backwards compatibility
        amrex::ParmParse pp_boussinesq_buoyancy("BoussinesqBuoyancy");
        amrex::ParmParse pp_abl("ABL");
        if (pp.contains("thermal_expansion_coefficient")) {
            pp.get("thermal_expansion_coefficient", m_constant_beta);
            if (pp_boussinesq_buoyancy.contains("thermal_expansion_coeff")) {
                amrex::Print()
                    << "WARNING: BoussinesqBuoyancy.thermal_expansion_coeff "
                       "option has been deprecated in favor of "
                       "transport.thermal_expansion_coefficient. Ignoring the "
                       "BoussinesqBuoyancy option in favor of the transport "
                       "option."
                    << std::endl;
            }
        } else if (pp_boussinesq_buoyancy.contains("thermal_expansion_coeff")) {
            amrex::Print()
                << "WARNING: BoussinesqBuoyancy.thermal_expansion_coeff option "
                   "has been deprecated in favor of "
                   "transport.thermal_expansion_coefficient. Please replace "
                   "this option."
                << std::endl;
            pp_boussinesq_buoyancy.get(
                "thermal_expansion_coeff", m_constant_beta);
        }

        if (pp.contains("reference_temperature")) {
            pp.get("reference_temperature", m_reference_temperature);
            if (pp_boussinesq_buoyancy.contains("reference_temperature")) {
                amrex::Print()
                    << "WARNING: BoussinesqBuoyancy.reference_temperature "
                       "option has been deprecated in favor of "
                       "transport.reference_temperature. Ignoring the "
                       "BoussinesqBuoyancy option in favor of the transport "
                       "option."
                    << std::endl;
            } else if (pp_abl.contains("reference_temperature")) {
                amrex::Print()
                    << "WARNING: ABL.reference_temperature "
                       "option has been deprecated in favor of "
                       "transport.reference_temperature. Ignoring the "
                       "ABL option in favor of the transport "
                       "option."
                    << std::endl;
            }
        } else if (pp_boussinesq_buoyancy.contains("reference_temperature")) {
            amrex::Print()
                << "WARNING: BoussinesqBuoyancy.reference_temperature option "
                   "has been deprecated in favor of "
                   "transport.reference_temperature. Please replace "
                   "this option."
                << std::endl;
            pp_boussinesq_buoyancy.get(
                "reference_temperature", m_reference_temperature);
        } else if (pp_abl.contains("reference_temperature")) {
            amrex::Print() << "WARNING: ABL.reference_temperature option "
                              "has been deprecated in favor of "
                              "transport.reference_temperature. Please replace "
                              "this option."
                           << std::endl;
            pp_abl.get("reference_temperature", m_reference_temperature);
        }
        if ((m_reference_temperature > 0.0) && (m_constant_beta > 0.0)) {
            amrex::Print()
                << "WARNING: specified both thermal_expansion_coefficient and "
                   "reference_temperature. Using reference_temperature to "
                   "define the thermal_expansion_coefficient. Remove "
                   "reference_temperature from the inputs if desiring "
                   "thermal_expansion_coefficient."
                << std::endl;
        }
    }

    ~TwoPhaseTransport() override = default;

    inline amrex::Real laminar_prandtl1() const { return m_Pr1; }
    inline amrex::Real laminar_prandtl2() const { return m_Pr2; }

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

            const auto& vof = m_repo.get_field("vof");

            for (int lev = 0; lev < m_repo.num_active_levels(); ++lev) {
                for (amrex::MFIter mfi((*mu)(lev)); mfi.isValid(); ++mfi) {
                    const auto& vbx = mfi.growntilebox();
                    const auto& volfrac = vof(lev).const_array(mfi);
                    const auto& visc = (*mu)(lev).array(mfi);
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

            const auto& levelset = m_repo.get_field("levelset");
            const auto& geom = m_repo.mesh().Geom();

            for (int lev = 0; lev < m_repo.num_active_levels(); ++lev) {
                for (amrex::MFIter mfi((*mu)(lev)); mfi.isValid(); ++mfi) {
                    const auto& vbx = mfi.growntilebox();
                    const auto& dx = geom[lev].CellSizeArray();
                    const auto& visc = (*mu)(lev).array(mfi);
                    const auto& phi = levelset(lev).const_array(mfi);
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

    //! Return the thermal diffusivity field (later divided by density, though)
    inline std::unique_ptr<ScratchField> alpha() override
    {
        // Select the interface capturing method
        auto alpha = m_repo.create_scratch_field(1, 1);

        if (m_ifacetype == InterfaceCapturingMethod::VOF) {

            const auto& vof = m_repo.get_field("vof");

            for (int lev = 0; lev < m_repo.num_active_levels(); ++lev) {
                for (amrex::MFIter mfi((*alpha)(lev)); mfi.isValid(); ++mfi) {
                    const auto& vbx = mfi.growntilebox();
                    const auto& volfrac = vof(lev).const_array(mfi);
                    const auto& thdiff = (*alpha)(lev).array(mfi);
                    const amrex::Real mu1 = m_mu1;
                    const amrex::Real mu2 = m_mu2;
                    const amrex::Real Pr1 = m_Pr1;
                    const amrex::Real Pr2 = m_Pr2;
                    amrex::ParallelFor(
                        vbx,
                        [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                            thdiff(i, j, k) =
                                mu1 / Pr1 * volfrac(i, j, k) +
                                mu2 / Pr2 * (1.0 - volfrac(i, j, k));
                        });
                }
            }

        } else if (m_ifacetype == InterfaceCapturingMethod::LS) {

            const auto& levelset = m_repo.get_field("levelset");
            const auto& geom = m_repo.mesh().Geom();

            for (int lev = 0; lev < m_repo.num_active_levels(); ++lev) {
                for (amrex::MFIter mfi((*alpha)(lev)); mfi.isValid(); ++mfi) {
                    const auto& vbx = mfi.growntilebox();
                    const auto& dx = geom[lev].CellSizeArray();
                    const auto& visc = (*alpha)(lev).array(mfi);
                    const auto& phi = levelset(lev).const_array(mfi);
                    const amrex::Real eps =
                        std::cbrt(2. * dx[0] * dx[1] * dx[2]);
                    const amrex::Real mu1 = m_mu1;
                    const amrex::Real mu2 = m_mu2;
                    const amrex::Real Pr1 = m_Pr1;
                    const amrex::Real Pr2 = m_Pr2;
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
                            visc(i, j, k) =
                                mu1 / Pr1 * smooth_heaviside +
                                mu2 / Pr2 * (1.0 - smooth_heaviside);
                        });
                }
            }
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

    //! Return the thermal expansion coefficient
    inline std::unique_ptr<ScratchField> beta() const override
    {
        auto beta = m_repo.create_scratch_field(1, 1);

        const amrex::Real beta_val = (m_reference_temperature > 0.0)
                                         ? 1.0 / m_reference_temperature
                                         : m_constant_beta;

        if (m_ifacetype == InterfaceCapturingMethod::VOF) {
            const auto& vof = m_repo.get_field("vof");
            for (int lev = 0; lev < m_repo.num_active_levels(); ++lev) {
                const auto& beta_arrs = (*beta)(lev).arrays();
                const auto& vof_arrs = vof(lev).const_arrays();
                amrex::ParallelFor(
                    (*beta)(lev), beta->num_grow(),
                    [=] AMREX_GPU_DEVICE(
                        int nbx, int i, int j, int k) noexcept {
                        if (vof_arrs[nbx](i, j, k) > constants::TIGHT_TOL) {
                            beta_arrs[nbx](i, j, k) = 0.0;
                        } else {
                            beta_arrs[nbx](i, j, k) = beta_val;
                        }
                    });
            }
        } else if (m_ifacetype == InterfaceCapturingMethod::LS) {
            const auto& levelset = m_repo.get_field("levelset");
            for (int lev = 0; lev < m_repo.num_active_levels(); ++lev) {
                const auto& dx = m_repo.mesh().Geom(lev).CellSizeArray();
                const auto& beta_arrs = (*beta)(lev).arrays();
                const auto& phi_arrs = levelset(lev).const_arrays();
                const amrex::Real eps = std::cbrt(2. * dx[0] * dx[1] * dx[2]);
                amrex::ParallelFor(
                    (*beta)(lev), beta->num_grow(),
                    [=] AMREX_GPU_DEVICE(
                        int nbx, int i, int j, int k) noexcept {
                        amrex::Real smooth_heaviside;
                        if (phi_arrs[nbx](i, j, k) > eps) {
                            smooth_heaviside = 1.0;
                        } else if (phi_arrs[nbx](i, j, k) < -eps) {
                            smooth_heaviside = 0.;
                        } else {
                            smooth_heaviside =
                                0.5 *
                                (1.0 + phi_arrs[nbx](i, j, k) / eps +
                                 1.0 / M_PI *
                                     std::sin(
                                         phi_arrs[nbx](i, j, k) * M_PI / eps));
                        }
                        beta_arrs[nbx](i, j, k) =
                            beta_val * (1 - smooth_heaviside);
                    });
            }
        }
        return beta;
    }

    inline amrex::Real reference_temperature() const override
    {
        return m_reference_temperature;
    }

    //! Return the reference temperature
    inline std::unique_ptr<ScratchField> ref_theta() const override
    {
        amrex::Abort("ref temp not implemented");
        auto ref_theta = m_repo.create_scratch_field(1, 1);
        return ref_theta;
    }

private:
    //! Reference to the CFD sim
    const CFDSim& m_sim;

    //! Reference to the field repository (for creating scratch fields)
    FieldRepo& m_repo;

    //! Interface capturing method variable
    InterfaceCapturingMethod m_ifacetype;

    //! Phase 1 (liquid) dynamic molecular viscosity
    amrex::Real m_mu1{1.0e-3};

    //! Phase 2 (gas) dynamic molecular viscosity
    amrex::Real m_mu2{1.0e-5};

    //! Phase 1 (liquid) Prandtl number
    amrex::Real m_Pr1{7.2};

    //! Phase 2 (gas) Prandtl number
    amrex::Real m_Pr2{0.7};

    //! Turbulent Prandtl number
    amrex::Real m_Prt{1.0};

    //! Constant thermal expansion coefficient
    amrex::Real m_constant_beta{0.0};

    //! Reference temperature
    amrex::Real m_reference_temperature{-1.0};
};

} // namespace amr_wind::transport

#endif /* TWOPHASETRANSPORT_H */
