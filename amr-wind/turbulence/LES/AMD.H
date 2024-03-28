#ifndef AMD_H
#define AMD_H

#include <string>
#include "amr-wind/turbulence/TurbModelBase.H"
#include "amr-wind/core/FieldRepo.H"
#include "amr-wind/fvm/stencils.H"

namespace amr_wind::turbulence {
/** AMD LES Model
 *  \ingroup turb_model
 */
template <typename Transport>
class AMD : public TurbModelBase<Transport>
{
public:
    static std::string identifier() { return "AMD-" + Transport::identifier(); }

    explicit AMD(CFDSim& sim);

    //! Model name for debugging purposes
    std::string model_name() const override { return "AMD"; }

    //! Update the turbulent viscosity field
    void update_turbulent_viscosity(
        const FieldState fstate, const DiffusionType /*unused*/) override;

    //! Update the effective thermal diffusivity field
    void update_alphaeff(Field& alphaeff) override;

    //! Return model coefficients dictionary
    TurbulenceModel::CoeffsDictType model_coeffs() const override;

    //! Parse turbulence model coefficients for this model
    void parse_model_coeffs() override;

    //! No post advance work for this model
    void post_advance_work() override{};

private:
    //! Poincare coefficient (default value set for 2nd order AMR-wind
    //! discretization)
    amrex::Real m_C{0.333333333333333};

    //! Reference temperature (Kelvin)
    amrex::Real m_ref_theta{300.0};

    const Field& m_vel;
    const Field& m_temperature;
    const Field& m_rho;
    amrex::Vector<amrex::Real> m_gravity{{0.0, 0.0, -9.81}};
};

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real amd_muvel(
    int i,
    int j,
    int k,
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dx, // Grid spacing
    amrex::Real beta, // Thermal expansion coefficient
    amrex::Real C,    // Poincare const
    amrex::Array4<amrex::Real const> const& gradVel,
    amrex::Array4<amrex::Real const> const& gradT) noexcept
{
    amrex::Real num_shear = 0;
    amrex::Real num_buoy = 0;
    amrex::Real denom = 0;
    for (int ii = 0; ii < AMREX_SPACEDIM; ++ii) {
        // num_buoy = dw/dx_i dT/dx_i
        num_buoy = num_buoy + gradVel(i, j, k, 2 * AMREX_SPACEDIM + ii) *
                                  gradT(i, j, k, ii) * dx[ii] * dx[ii];
        for (int jj = 0; jj < AMREX_SPACEDIM; ++jj) {
            denom = denom + gradVel(i, j, k, ii * AMREX_SPACEDIM + jj) *
                                gradVel(i, j, k, ii * AMREX_SPACEDIM + jj);
            amrex::Real sij =
                0.5 * (gradVel(i, j, k, ii * AMREX_SPACEDIM + jj) +
                       gradVel(i, j, k, jj * AMREX_SPACEDIM + ii));
            for (int kk = 0; kk < AMREX_SPACEDIM; ++kk) {
                amrex::Real dkui = gradVel(i, j, k, ii * AMREX_SPACEDIM + kk);
                amrex::Real dkuj = gradVel(i, j, k, jj * AMREX_SPACEDIM + kk);
                num_shear = num_shear + dkui * dkuj * dx[kk] * dx[kk] * sij;
            }
        }
    }
    denom = std::max(1e-15, denom);
    num_shear = -C * num_shear;
    num_buoy = C * beta * num_buoy;
    return std::max(1e-15, (num_shear + num_buoy) / denom);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real amd_thermal_diff(
    int i,
    int j,
    int k,
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dx, // Grid spacing
    amrex::Real C,    // Poincare const
    amrex::Array4<amrex::Real const> const& gradVel,
    amrex::Array4<amrex::Real const> const& gradT) noexcept
{
    amrex::Real num = 0;
    amrex::Real denom = 0;
    for (int ii = 0; ii < AMREX_SPACEDIM; ++ii) {
        amrex::Real diT = gradT(i, j, k, ii);
        denom = denom + diT * diT;
        for (int kk = 0; kk < AMREX_SPACEDIM; ++kk) {
            amrex::Real dkui = gradVel(i, j, k, ii * AMREX_SPACEDIM + kk);
            amrex::Real dkT = gradT(i, j, k, kk);
            num = num + dkui * diT * dkT * dx[kk] * dx[kk];
        }
    }

    denom = std::max(1e-15, denom);
    num = -C * num;
    return std::max(1e-15, num / denom);
}

} // namespace amr_wind::turbulence

#endif /* AMD_H */
