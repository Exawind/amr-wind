#ifndef AMD_H
#define AMD_H

#include <AMReX_GpuContainers.H>
#include <string>
#include "amr-wind/turbulence/TurbModelBase.H"
#include "amr-wind/utilities/FieldPlaneAveraging.H"
#include "amr-wind/core/FieldRepo.H"
#include "amr-wind/fvm/stencils.H"
#include "amr-wind/utilities/DirectionSelector.H"
#include "amr-wind/utilities/linear_interpolation.H"

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

    template <typename IndexSelector>
    void update_turbulent_viscosity(
        const FieldState fstate, const IndexSelector& idxOp);

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

    //! Wall-normal direction axis
    int m_normal_dir{2};

    const Field& m_vel;
    const Field& m_temperature;
    const Field& m_rho;
    FieldPlaneAveraging m_pa_temp;
    amrex::Vector<amrex::Real> m_gravity{0.0, 0.0, -9.81};
};

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real amd_muvel(
    int i,
    int j,
    int k,
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& dx, // Grid spacing
    const amrex::Real beta, // Thermal expansion coefficient
    const amrex::Real C,    // Poincare const
    const amrex::Array4<amrex::Real const>& gradVel,
    const amrex::Array4<amrex::Real const>& gradT,
    const amrex::Real* gradTbar_coord_begin,
    const amrex::Real* gradTbar_coord_end,
    const amrex::Real* gradTbar,
    const int normal_dir,
    const amrex::Real nlo) noexcept
{
    amrex::Real num_shear = 0;
    amrex::Real num_buoy = 0;
    amrex::Real denom = 0;
    const int n_ind = amr_wind::direction_selector(i, j, k, normal_dir);
    const amrex::Real h = nlo + (n_ind + 0.5) * dx[normal_dir];
    const amrex::Real gradTbar_h = amr_wind::interp::linear(
        gradTbar_coord_begin, gradTbar_coord_end, gradTbar, h);
    for (int ii = 0; ii < AMREX_SPACEDIM; ++ii) {
        // This should operate only on the wall normal velocity
        num_buoy +=
            gradVel(i, j, k, normal_dir * AMREX_SPACEDIM + ii) *
            (gradT(i, j, k, ii) - ((ii == normal_dir) ? gradTbar_h : 0.0)) *
            dx[ii] * dx[ii];

        for (int jj = 0; jj < AMREX_SPACEDIM; ++jj) {
            const amrex::Real diuj = gradVel(i, j, k, ii * AMREX_SPACEDIM + jj);
            const amrex::Real djui = gradVel(i, j, k, jj * AMREX_SPACEDIM + ii);
            denom += diuj * diuj;
            const amrex::Real sij = 0.5 * (diuj + djui);
            for (int kk = 0; kk < AMREX_SPACEDIM; ++kk) {
                const amrex::Real dkui =
                    gradVel(i, j, k, ii * AMREX_SPACEDIM + kk);
                const amrex::Real dkuj =
                    gradVel(i, j, k, jj * AMREX_SPACEDIM + kk);
                num_shear += dkui * dkuj * dx[kk] * dx[kk] * sij;
            }
        }
    }
    denom = std::max(1e-15, denom);
    num_shear *= -C;
    num_buoy *= C * beta;
    return std::max(1e-15, (num_shear + num_buoy) / denom);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real amd_thermal_diff(
    int i,
    int j,
    int k,
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& dx, // Grid spacing
    const amrex::Real C,                                    // Poincare const
    const amrex::Array4<amrex::Real const>& gradVel,
    const amrex::Array4<amrex::Real const>& gradT) noexcept
{
    amrex::Real num = 0;
    amrex::Real denom = 0;
    for (int ii = 0; ii < AMREX_SPACEDIM; ++ii) {
        const amrex::Real diT = gradT(i, j, k, ii);
        denom += diT * diT;
        for (int kk = 0; kk < AMREX_SPACEDIM; ++kk) {
            const amrex::Real dkui = gradVel(i, j, k, ii * AMREX_SPACEDIM + kk);
            const amrex::Real dkT = gradT(i, j, k, kk);
            num += dkui * diT * dkT * dx[kk] * dx[kk];
        }
    }

    denom = std::max(1e-15, denom);
    num *= -C;
    return std::max(1e-15, num / denom);
}

} // namespace amr_wind::turbulence

#endif /* AMD_H */
