#ifndef CONSTTRANSPORT_H
#define CONSTTRANSPORT_H

#include "amr-wind/transport_models/TransportModel.H"
#include "amr-wind/utilities/constants.H"
#include "AMReX_ParmParse.H"

namespace amr_wind::transport {

/** Fluid transport model with constant properties
 *  \ingroup transport
 */
class ConstTransport : public TransportModel::Register<ConstTransport>
{
public:
    static constexpr bool constant_properties = true;

    static std::string identifier() { return "ConstTransport"; }

    explicit ConstTransport(const CFDSim& sim)
        : m_repo(sim.repo()), m_phy_mgr(sim.physics_manager())
    {
        amrex::ParmParse pp("transport");
        pp.query("viscosity", m_mu);
        pp.query("laminar_prandtl", m_Pr);
        pp.query("turbulent_prandtl", m_Prt);

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

    ~ConstTransport() override = default;

    inline amrex::Real viscosity() const { return m_mu; }

    inline amrex::Real thermal_diffusivity() const { return m_mu / m_Pr; }

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
        auto mu = m_repo.create_scratch_field(1, 1);
        for (int lev = 0; lev < m_repo.num_active_levels(); ++lev) {
            (*mu)(lev).setVal(m_mu);
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

    //! Return the thermal expansion coefficient
    inline std::unique_ptr<ScratchField> beta() const override
    {
        auto beta = m_repo.create_scratch_field(1, 1);

        const amrex::Real beta_val = (m_reference_temperature > 0.0)
                                         ? 1.0 / m_reference_temperature
                                         : m_constant_beta;

        for (int lev = 0; lev < m_repo.num_active_levels(); ++lev) {
            (*beta)(lev).setVal(beta_val);
        }

        if (m_repo.field_exists("vof")) {
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
                        }
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
        if (m_reference_temperature < 0.0) {
            amrex::Abort("Reference temperature was not set");
        }

        auto ref_theta = m_repo.create_scratch_field(1, 1);
        for (int lev = 0; lev < m_repo.num_active_levels(); ++lev) {
            (*ref_theta)(lev).setVal(m_reference_temperature);
        }
        return ref_theta;
    }

private:
    //! Reference to the field repository (for creating scratch fields)
    FieldRepo& m_repo;

    //! Reference to the physics manager
    const PhysicsMgr& m_phy_mgr;

    //! (Laminar) dynamic viscosity
    amrex::Real m_mu{1.0e-5};

    //! Prandtl number
    amrex::Real m_Pr{1.0};

    //! Turbulent Prandtl number
    amrex::Real m_Prt{1.0};

    //! Constant thermal expansion coefficient
    amrex::Real m_constant_beta{0.0};

    //! Reference temperature
    amrex::Real m_reference_temperature{-1.0};
};

} // namespace amr_wind::transport

#endif /* CONSTTRANSPORT_H */
