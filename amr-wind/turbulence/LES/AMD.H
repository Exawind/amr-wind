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

template <typename Stencil>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real amd_muvel(
    int i,
    int j,
    int k,
    amrex::Real dx,   // Grid spacing in x direction
    amrex::Real dy,   // Grid spacing in x direction
    amrex::Real dz,   // Grid spacing in x direction
    amrex::Real beta, // Thermal expansion coefficient
    amrex::Real C,    // Poincare const
    amrex::Array4<amrex::Real const> const& vel,
    amrex::Array4<amrex::Real const> const& theta) noexcept
{
    const amrex::Real idx = 1.0 / dx;
    const amrex::Real idy = 1.0 / dy;
    const amrex::Real idz = 1.0 / dz;

    amrex::Real cp1, c, cm1;

    cp1 = Stencil::c00;
    c = Stencil::c01;
    cm1 = Stencil::c02;

    const amrex::Real ux = (cp1 * vel(i + 1, j, k, 0) + c * vel(i, j, k, 0) +
                            cm1 * vel(i - 1, j, k, 0)) *
                           idx;
    const amrex::Real vx = (cp1 * vel(i + 1, j, k, 1) + c * vel(i, j, k, 1) +
                            cm1 * vel(i - 1, j, k, 1)) *
                           idx;
    const amrex::Real wx = (cp1 * vel(i + 1, j, k, 2) + c * vel(i, j, k, 2) +
                            cm1 * vel(i - 1, j, k, 2)) *
                           idx;
    const amrex::Real tx = (cp1 * theta(i + 1, j, k) + c * theta(i, j, k) +
                            cm1 * theta(i - 1, j, k)) *
                           idx;

    cp1 = Stencil::c10;
    c = Stencil::c11;
    cm1 = Stencil::c12;

    const amrex::Real uy = (cp1 * vel(i, j + 1, k, 0) + c * vel(i, j, k, 0) +
                            cm1 * vel(i, j - 1, k, 0)) *
                           idy;
    const amrex::Real vy = (cp1 * vel(i, j + 1, k, 1) + c * vel(i, j, k, 1) +
                            cm1 * vel(i, j - 1, k, 1)) *
                           idy;
    const amrex::Real wy = (cp1 * vel(i, j + 1, k, 2) + c * vel(i, j, k, 2) +
                            cm1 * vel(i, j - 1, k, 2)) *
                           idy;
    const amrex::Real ty = (cp1 * theta(i, j + 1, k) + c * theta(i, j, k) +
                            cm1 * theta(i, j - 1, k)) *
                           idy;

    cp1 = Stencil::c20;
    c = Stencil::c21;
    cm1 = Stencil::c22;

    const amrex::Real uz = (cp1 * vel(i, j, k + 1, 0) + c * vel(i, j, k, 0) +
                            cm1 * vel(i, j, k - 1, 0)) *
                           idz;
    const amrex::Real vz = (cp1 * vel(i, j, k + 1, 1) + c * vel(i, j, k, 1) +
                            cm1 * vel(i, j, k - 1, 1)) *
                           idz;
    const amrex::Real wz = (cp1 * vel(i, j, k + 1, 2) + c * vel(i, j, k, 2) +
                            cm1 * vel(i, j, k - 1, 2)) *
                           idz;
    const amrex::Real tz = (cp1 * theta(i, j, k + 1) + c * theta(i, j, k) +
                            cm1 * theta(i, j, k - 1)) *
                           idz;

    const amrex::Real dx2 = dx * dx;
    const amrex::Real dy2 = dy * dy;
    const amrex::Real dz2 = dz * dz;
    const amrex::Real num_shear =
        -2.0 * C *
        (ux * (ux * ux * dx2 + uy * uy * dy2 + uz * uz * dz2) +
         vy * (vx * vx * dx2 + vy * vy * dy2 + vz * vz * dz2) +
         wz * (wx * wx * dx2 + wy * wy * dy2 + wz * wz * dz2) +
         (uy + vx) * (ux * vx * dx2 + uy * vy * dy2 + uz * vz * dz2) +
         (uz + wx) * (ux * wx * dx2 + uy * wy * dy2 + uz * wz * dz2) +
         (vz + wy) * (vx * wx * dx2 + vy * wy * dy2 + vz * wz * dz2));
    amrex::Real num_buoy =
        C * beta * (dx2 * wx * tx + dy2 * wy * ty + dz2 * wz * tz);
    amrex::Real denom = (ux * ux + uy * uy + uz * uz + vx * vx + vy * vy +
                         vz * vz + wx * wx + wy * wy + wz * wz) +
                        1e-15;
    return std::max(1e-15, (num_shear + num_buoy) / denom);
}

template <typename Stencil>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real amd_thermal_diff(
    int i,
    int j,
    int k,
    amrex::Real dx, // Grid spacing in x direction
    amrex::Real dy, // Grid spacing in x direction
    amrex::Real dz, // Grid spacing in x direction
    amrex::Real C,  // Poincare const
    amrex::Array4<amrex::Real const> const& vel,
    amrex::Array4<amrex::Real const> const& theta) noexcept
{
    const amrex::Real idx = 1.0 / dx;
    const amrex::Real idy = 1.0 / dy;
    const amrex::Real idz = 1.0 / dz;

    amrex::Real cp1, c, cm1;

    cp1 = Stencil::c00;
    c = Stencil::c01;
    cm1 = Stencil::c02;

    const amrex::Real ux = (cp1 * vel(i + 1, j, k, 0) + c * vel(i, j, k, 0) +
                            cm1 * vel(i - 1, j, k, 0)) *
                           idx;
    const amrex::Real vx = (cp1 * vel(i + 1, j, k, 1) + c * vel(i, j, k, 1) +
                            cm1 * vel(i - 1, j, k, 1)) *
                           idx;
    const amrex::Real wx = (cp1 * vel(i + 1, j, k, 2) + c * vel(i, j, k, 2) +
                            cm1 * vel(i - 1, j, k, 2)) *
                           idx;
    const amrex::Real tx = (cp1 * theta(i + 1, j, k) + c * theta(i, j, k) +
                            cm1 * theta(i - 1, j, k)) *
                           idx;

    cp1 = Stencil::c10;
    c = Stencil::c11;
    cm1 = Stencil::c12;

    const amrex::Real uy = (cp1 * vel(i, j + 1, k, 0) + c * vel(i, j, k, 0) +
                            cm1 * vel(i, j - 1, k, 0)) *
                           idy;
    const amrex::Real vy = (cp1 * vel(i, j + 1, k, 1) + c * vel(i, j, k, 1) +
                            cm1 * vel(i, j - 1, k, 1)) *
                           idy;
    const amrex::Real wy = (cp1 * vel(i, j + 1, k, 2) + c * vel(i, j, k, 2) +
                            cm1 * vel(i, j - 1, k, 2)) *
                           idy;
    const amrex::Real ty = (cp1 * theta(i, j + 1, k) + c * theta(i, j, k) +
                            cm1 * theta(i, j - 1, k)) *
                           idy;

    cp1 = Stencil::c20;
    c = Stencil::c21;
    cm1 = Stencil::c22;

    const amrex::Real uz = (cp1 * vel(i, j, k + 1, 0) + c * vel(i, j, k, 0) +
                            cm1 * vel(i, j, k - 1, 0)) *
                           idz;
    const amrex::Real vz = (cp1 * vel(i, j, k + 1, 1) + c * vel(i, j, k, 1) +
                            cm1 * vel(i, j, k - 1, 1)) *
                           idz;
    const amrex::Real wz = (cp1 * vel(i, j, k + 1, 2) + c * vel(i, j, k, 2) +
                            cm1 * vel(i, j, k - 1, 2)) *
                           idz;
    const amrex::Real tz = (cp1 * theta(i, j, k + 1) + c * theta(i, j, k) +
                            cm1 * theta(i, j, k - 1)) *
                           idz;

    const amrex::Real dx2 = dx * dx;
    const amrex::Real dy2 = dy * dy;
    const amrex::Real dz2 = dz * dz;
    const amrex::Real num =
        -C * ((dx2 * ux * tx + dy2 * uy * ty + dz2 * uz * tz) * tx +
              (dx2 * vx * tx + dy2 * vy * ty + dz2 * vz * tz) * ty +
              (dx2 * wx * tx + dy2 * wy * ty + dz2 * wz * tz) * tz);

    amrex::Real denom = tx * tx + ty * ty + tz * tz + 1e-15;
    return std::max(1e-15, num / denom);
}

} // namespace amr_wind::turbulence

#endif /* AMD_H */
