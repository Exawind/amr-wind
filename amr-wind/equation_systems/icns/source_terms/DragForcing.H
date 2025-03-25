#ifndef DRAGFORCING_H
#define DRAGFORCING_H

#include "amr-wind/equation_systems/icns/MomentumSource.H"
#include "amr-wind/core/SimTime.H"
#include "amr-wind/CFDSim.H"

namespace amr_wind::pde::icns {

/** Adds the forcing term to include the presence of immersed boundary
 *
 *  \ingroup icns_src
 *
 *
 */
class DragForcing : public MomentumSource::Register<DragForcing>
{
public:
    static std::string identifier() { return "DragForcing"; }

    explicit DragForcing(const CFDSim& sim);

    ~DragForcing() override;

    void operator()(
        const int lev,
        const amrex::MFIter& mfi,
        const amrex::Box& bx,
        const FieldState fstate,
        const amrex::Array4<amrex::Real>& src_term) const override;

private:
    const SimTime& m_time;
    const CFDSim& m_sim;
    const amrex::AmrCore& m_mesh;
    const Field& m_velocity;
    const Field* m_target_vel{nullptr};
    const Field* m_target_levelset{nullptr};
    amrex::Gpu::DeviceVector<amrex::Real> m_device_vel_ht;
    amrex::Gpu::DeviceVector<amrex::Real> m_device_vel_vals;
    amrex::Real m_drag_coefficient{10.0};
    amrex::Real m_sponge_strength{1.0};
    amrex::Real m_sponge_density{1.0};
    amrex::Real m_sponge_distance_west{-1000};
    amrex::Real m_sponge_distance_east{1000};
    amrex::Real m_sponge_distance_south{-1000};
    amrex::Real m_sponge_distance_north{1000};
    int m_sponge_west{0};
    int m_sponge_east{1};
    int m_sponge_south{0};
    int m_sponge_north{1};
    bool m_is_laminar{false};
    bool m_terrain_is_waves{false};
    std::string m_wall_het_model{"none"};
    amrex::Real m_mol_length{-1e30};
    amrex::Real m_kappa{0.41};
    amrex::Real m_gamma_m{5.0};
    amrex::Real m_beta_m{16.0};

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
};

} // namespace amr_wind::pde::icns

#endif
