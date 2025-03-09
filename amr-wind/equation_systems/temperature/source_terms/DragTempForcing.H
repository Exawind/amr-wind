#ifndef DRAGTEMPFORCING_H
#define DRAGTEMPFORCING_H

#include "amr-wind/equation_systems/temperature/TemperatureSource.H"
#include "amr-wind/core/SimTime.H"
#include "amr-wind/CFDSim.H"

namespace amr_wind::pde::temperature {

class DragTempForcing : public TemperatureSource::Register<DragTempForcing>
{
public:
    static std::string identifier() { return "DragTempForcing"; }

    explicit DragTempForcing(const CFDSim& sim);

    ~DragTempForcing() override;

    void operator()(
        const int lev,
        const amrex::MFIter& mfi,
        const amrex::Box& bx,
        const FieldState /*fstate*/,
        const amrex::Array4<amrex::Real>& src_term) const override;

private:
    const SimTime& m_time;
    const CFDSim& m_sim;
    const amrex::AmrCore& m_mesh;
    const Field& m_velocity;
    const Field& m_temperature;
    amrex::Real m_drag_coefficient{1.0};

    std::string m_wall_het_model{"none"};
    amrex::Real m_mol_length{-1e30};
    amrex::Real m_kappa{0.41};
    amrex::Real m_gamma_m{5.0};
    amrex::Real m_beta_m{16.0};
    amrex::Real m_gamma_h{5.0};
    amrex::Real m_beta_h{16.0};
    amrex::Real m_soil_temperature{300};
    amrex::Vector<amrex::Real> m_gravity{0.0, 0.0, -9.81};

    static amrex::Real stability(
        const amrex::Real z,
        const amrex::Real mol_length,
        const amrex::Real gamma_m,
        const amrex::Real beta_m)
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

    static amrex::Real thermal_stability(
        const amrex::Real z,
        const amrex::Real mol_length,
        const amrex::Real gamma_h,
        const amrex::Real beta_h)
    {
        const amrex::Real zeta = z / mol_length;
        amrex::Real psi_h = 0;
        if (zeta > 0) {
            psi_h = -gamma_h * zeta;
        } else {
            amrex::Real x = std::sqrt(1 - beta_h * zeta);
            psi_h = 2.0 * std::log(0.5 * (1 + x));
        }
        return psi_h;
    }
};

} // namespace amr_wind::pde::temperature
#endif
