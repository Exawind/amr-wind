#ifndef AMDNOTHERM_H
#define AMDNOTHERM_H

#include <AMReX_Config.H>
#include <AMReX_AmrCore.H>
#include <string>
#include "amr-wind/incflo_enums.H"
#include "amr-wind/turbulence/TurbModelBase.H"
#include "amr-wind/core/FieldRepo.H"
#include "amr-wind/fvm/stencils.H"

namespace amr_wind::turbulence {
/** AMD LES Model
 *  \ingroup turb_model
 */
template <typename Transport>
class AMDNoTherm : public TurbModelBase<Transport>
{
public:
    static std::string identifier()
    {
        return "AMDNoTherm-" + Transport::identifier();
    }

    explicit AMDNoTherm(CFDSim& sim);

    //! Model name for debugging purposes
    std::string model_name() const override { return "AMDNoTherm"; }

    //! Update the turbulent viscosity field
    void update_turbulent_viscosity(
        const FieldState fstate, const DiffusionType /*unused*/) override;

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
    const Field& m_vel;
    const Field& m_rho;
};

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real amd_base_muvel(
    int i,
    int j,
    int k,
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dx, // Grid spacing
    amrex::Real C,                                         // Poincare const
    amrex::Array4<amrex::Real const> const& gradVel) noexcept
{

    amrex::Real num_shear = 0;
    amrex::Real denom = 0;
    for (int ii = 0; ii < AMREX_SPACEDIM; ++ii) {
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

    return std::max(1e-15, num_shear / denom);
}

} // namespace amr_wind::turbulence

#endif /* AMDNoTherm_H */
