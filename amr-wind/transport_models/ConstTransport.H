#ifndef CONSTTRANSPORT_H
#define CONSTTRANSPORT_H

#include "amr-wind/transport_models/TransportModel.H"
#include "AMReX_ParmParse.H"

namespace amr_wind::transport {

/** Fluid transport model with constant properties
 *  \ingroup transport
 */
class ConstTransport : public TransportModel
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
        if (pp.contains("thermal_expansion_coefficient")) {
            pp.get("thermal_expansion_coefficient", m_beta);
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
            pp_boussinesq_buoyancy.get("thermal_expansion_coeff", m_beta);
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

    //! Return the thermal diffusivity field
    inline std::unique_ptr<ScratchField> beta() override
    {
        auto beta = m_repo.create_scratch_field(1, 1);

        amrex::Real beta_const = m_beta;
        if (m_phy_mgr.contains("ABL")) {
            const int normal_dir = 2;
            amrex::Vector<amrex::Real> gravity{0.0, 0.0, -9.81};
            amrex::Real ref_theta = 0.0;
            {
                amrex::ParmParse pp("ABL");
                pp.get("reference_temperature", ref_theta);
            }
            {
                amrex::ParmParse pp("incflo");
                pp.queryarr("gravity", gravity);
            }
            beta_const = -gravity[normal_dir] / ref_theta;
        }

        for (int lev = 0; lev < m_repo.num_active_levels(); ++lev) {
            (*beta)(lev).setVal(beta_const);
        }

        return beta;
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

    //! Thermal expansion coefficient
    amrex::Real m_beta{0.0};
};

} // namespace amr_wind::transport

#endif /* CONSTTRANSPORT_H */
