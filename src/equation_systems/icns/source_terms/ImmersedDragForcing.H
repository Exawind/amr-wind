#ifndef IMMERSEDDRAGFORCING_H
#define IMMERSEDDRAGFORCING_H

#include "src/equation_systems/icns/MomentumSource.H"
#include "src/core/SimTime.H"
#include "src/CFDSim.H"
#include "src/utilities/constants.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace kynema_sgf::pde::icns {

/** Momentum forcing for terrain represented by the ImmersedTerrain physics
 *
 *  \ingroup icns_src
 *
 *  Consumes ``terrain_fraction``, ``terrain_mask``, ``terrain_surface`` and
 *  ``terrain_roughness`` declared by ImmersedTerrain. Two contributions:
 *
 *  1. Immersed drag in every cell with a non-zero terrain fraction beta,
 *     relaxing the velocity toward zero at the rate
 *     C = beta * drag_coefficient / dz. The relaxation is integrated exactly
 *     over the time step, C_eff = (1 - exp(-C dt)) / dt, so it is stable and
 *     monotone for any C and dt with no ad hoc limiters.
 *
 *  2. Wall model (MOST log law) in surface cells (mask == 2), weighted by
 *     (1 - beta), in the turbulent pathway only. Three ways of measuring the
 *     wall distance are selectable with ``wall_model``:
 *       - ``cell_offset``    : d1 = 0.5 dx_f, d2 = 1.5 dx_f on each of the
 *                              six faces that touch a solid neighbour
 *       - ``terrain_height`` : as above, but the bottom face uses the true
 *                              height above the terrain, z - h
 *       - ``surface_normal`` : distance and velocity decomposition along the
 *                              terrain normal from the surface slopes; the
 *                              stress acts on u, v and w together. Side-wall
 *                              cells without terrain below fall back to the
 *                              face-based treatment.
 *     The wall model applies the stress divergence -u*^2 e_t / dx_n and a
 *     relaxation of the tangential velocity toward the log-law value with
 *     time scale ``bc_forcing_time_scale``: ``wall`` uses
 *     max(tau_f dt, d1/u*), ``time_step`` uses tau_f dt.
 *
 *  The pathway is selected from ``turbulence.model``: ``Laminar`` gives drag
 *  only, anything else adds the wall model. ``force_laminar`` overrides.
 */
class ImmersedDragForcing : public MomentumSource::Register<ImmersedDragForcing>
{
public:
    static std::string identifier() { return "ImmersedDragForcing"; }

    explicit ImmersedDragForcing(const CFDSim& sim);

    ~ImmersedDragForcing() override;

    void operator()(
        int lev, FieldState fstate, amrex::MultiFab& src_term) const override;

private:
    const SimTime& m_time;
    const CFDSim& m_sim;
    const amrex::AmrCore& m_mesh;
    const Field& m_velocity;

    //! Immersed drag relaxation rate is drag_coefficient / dz
    amrex::Real m_drag_coefficient{10.0_rt};
    //! Wall model relaxation time scale in units of dt (floor for ``wall``)
    amrex::Real m_forcing_time_factor{5.0_rt};
    //! ``wall`` or ``time_step``
    std::string m_bc_forcing_time_scale{"wall"};
    //! ``cell_offset``, ``terrain_height`` or ``surface_normal``
    std::string m_wall_model{"cell_offset"};
    //! Lower bound on roughness length
    amrex::Real m_min_z0{1.0e-4_rt};
    //! Neighbour fraction at or above which a face is treated as a wall
    amrex::Real m_solid_threshold{0.5_rt};
    //! Skip the wall model even with a turbulence model active
    bool m_force_laminar{false};
    bool m_is_laminar{false};

    //! Monin-Obukhov inputs shared with the ABL wall function
    std::string m_wall_het_model{"none"};
    amrex::Real m_monin_obukhov_length{constants::LARGE_NUM};
    amrex::Real m_kappa{0.41_rt};
    amrex::Real m_gamma_m{5.0_rt};
    amrex::Real m_beta_m{16.0_rt};
};

} // namespace kynema_sgf::pde::icns

#endif
