#ifndef IMMERSEDDRAGTEMPFORCING_H
#define IMMERSEDDRAGTEMPFORCING_H

#include "src/equation_systems/temperature/TemperatureSource.H"
#include "src/core/SimTime.H"
#include "src/CFDSim.H"
#include "src/utilities/constants.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace kynema_sgf::pde::temperature {

/** Temperature forcing for terrain represented by the ImmersedTerrain physics
 *
 *  Counterpart of ImmersedDragForcing for the temperature equation.
 *
 *  1. Inside the terrain (fraction beta > 0) the temperature relaxes toward
 *     ``soil_temperature`` at the rate beta * drag_coefficient / dz, integrated
 *     exactly over the step. Following Munoz-Esparza et al. (2020) the rate
 *     does not depend on the local velocity.
 *
 *  2. In surface cells (mask == 2), weighted by (1 - beta) and only when the
 *     turbulence model is not ``Laminar``, a MOST heat-flux wall model is
 *     applied on the same wall patches as the momentum source (shared
 *     ImmersedWallModel.H; the wall_model, time scale and roughness inputs
 *     are read from the ``ImmersedDragForcing`` prefix so the two cannot
 *     drift apart). The friction velocity comes from the reference cell at d2;
 *     the temperature scale theta* depends on ``surface_condition``:
 *       - ``obukhov_length``      : theta* = theta u*^2 / (kappa g L) with the
 *                                   single ABL Obukhov length (current
 *                                   DragTempForcing behavior)
 *       - ``surface_temperature`` : theta* = kappa (theta_ref - theta_soil) /
 *                                   phi_h(d2), surface held at soil_temperature
 *       - ``heat_flux``           : theta* = -surface_heat_flux / u*
 *     The forcing is the surface flux divergence q_s / dx_n plus a relaxation
 *     of the cell temperature toward the log-law value at d1.
 */
class ImmersedDragTempForcing
    : public TemperatureSource::Register<ImmersedDragTempForcing>
{
public:
    static std::string identifier() { return "ImmersedDragTempForcing"; }

    explicit ImmersedDragTempForcing(const CFDSim& sim);

    ~ImmersedDragTempForcing() override;

    void operator()(
        int lev, FieldState fstate, amrex::MultiFab& src_term) const override;

private:
    const SimTime& m_time;
    const CFDSim& m_sim;
    const amrex::AmrCore& m_mesh;
    const Field& m_velocity;
    const Field& m_temperature;

    //! Interior relaxation rate is drag_coefficient / dz
    amrex::Real m_drag_coefficient{10.0_rt};
    //! Temperature inside the terrain, and the surface temperature for
    //! ``surface_condition = surface_temperature``
    amrex::Real m_soil_temperature{300.0_rt};
    //! ``obukhov_length``, ``surface_temperature`` or ``heat_flux``
    std::string m_surface_condition{"obukhov_length"};
    //! Kinematic surface heat flux for ``surface_condition = heat_flux``
    amrex::Real m_surface_heat_flux{0.0_rt};
    //! Skip the wall model even with a turbulence model active
    bool m_force_laminar{false};
    bool m_is_laminar{false};

    //! Shared with ImmersedDragForcing (read from its input prefix)
    std::string m_wall_model{"cell_offset"};
    std::string m_reference_distance{"nominal"};
    std::string m_drag_weight{"fraction"};
    std::string m_bc_forcing_time_scale{"wall"};
    amrex::Real m_forcing_time_factor{5.0_rt};
    amrex::Real m_min_z0{1.0e-4_rt};
    amrex::Real m_solid_threshold{0.5_rt};

    //! Monin-Obukhov inputs shared with the ABL wall function
    std::string m_wall_het_model{"none"};
    amrex::Real m_monin_obukhov_length{constants::LARGE_NUM};
    amrex::Real m_kappa{0.41_rt};
    amrex::Real m_gamma_m{5.0_rt};
    amrex::Real m_beta_m{16.0_rt};
    amrex::Real m_gamma_h{5.0_rt};
    amrex::Real m_beta_h{16.0_rt};
    amrex::Vector<amrex::Real> m_gravity{0.0_rt, 0.0_rt, -9.81_rt};
};

} // namespace kynema_sgf::pde::temperature

#endif
