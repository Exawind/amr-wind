#ifndef KRANSAXELL_H
#define KRANSAXELL_H

#include "src/transport_models/TransportModel.H"
#include "src/equation_systems/tke/TKESource.H"
#include "src/utilities/constants.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace kynema_sgf::pde::tke {

/** TKE source term based on Axell 2011 paper
 * Axell, L. B., & Liungman, O. (2001). A one-equation turbulence model for
 * geophysical applications: comparison with data and the k− ε model.
 * Environmental Fluid Mechanics, 1, 71-106.
 *  \ingroup tke_src turb_model we_abl
 *
 *  With terrain (``KLAxell.terrain_model``): the source is blanked inside the
 *  terrain, the TKE is relaxed toward the log-law value
 *  (u*^3/C_mu^3 + B)^(2/3) in the first fluid cells and damped inside the
 *  terrain. ``TerrainDrag`` (default) uses the binary blanking and the drag
 *  cells; ``ImmersedTerrain`` evaluates the same forcing in a separate kernel
 *  from the terrain fraction and the wall patches of the surface cells shared
 *  with ImmersedDragForcing, with the interior damping integrated exactly in
 *  time.
 */
class KransAxell : public TKESource::Register<KransAxell>
{
public:
    static std::string identifier() { return "KransAxell"; }

    explicit KransAxell(const CFDSim& /*sim*/);

    ~KransAxell() override;

    void operator()(
        int lev, FieldState fstate, amrex::MultiFab& src_term) const override;

private:
    //! Terrain forcing on one level with the ImmersedTerrain fields
    //! (``KLAxell.terrain_model = ImmersedTerrain``)
    void immersed_terrain_forcing(
        int lev, FieldState fstate, amrex::MultiFab& src_term) const;

    Field& m_turb_lscale;
    Field& m_shear_prod;
    Field& m_buoy_prod;
    Field& m_dissip;
    Field& m_tke;
    amrex::Real m_Cmu{0.556_rt};
    amrex::Real m_heat_flux{0.0_rt};
    amrex::Real m_z0{0.1_rt};
    amrex::Real m_kappa{0.41_rt};
    amrex::Real m_meso_start{600.0_rt};
    amrex::Real m_forcing_time_factor{5.0_rt};
    amrex::Vector<amrex::Real> m_gravity{0.0_rt, 0.0_rt, -9.81_rt};
    const SimTime& m_time;
    const CFDSim& m_sim;
    const amrex::AmrCore& m_mesh;
    const Field& m_velocity;
    std::string m_1d_rans;
    amrex::Vector<amrex::Real> m_wind_heights;
    amrex::Vector<amrex::Real> m_tke_values;
    amrex::Gpu::DeviceVector<amrex::Real> m_wind_heights_d;
    amrex::Gpu::DeviceVector<amrex::Real> m_tke_values_d;
    //! Transport model
    const transport::TransportModel& m_transport;

    //! Reference temperature
    std::unique_ptr<ScratchField> m_ref_theta;
    mutable std::unique_ptr<ScratchField> m_ref_theta_scratch;
    bool m_horizontal_sponge{false};
    amrex::Real m_sponge_strength{1.0_rt};
    amrex::Real m_sponge_distance_west{0.0_rt};
    amrex::Real m_sponge_distance_east{0.0_rt};
    amrex::Real m_sponge_distance_south{0.0_rt};
    amrex::Real m_sponge_distance_north{0.0_rt};
    bool m_sponge_west{false};
    bool m_sponge_east{false};
    bool m_sponge_south{false};
    bool m_sponge_north{false};
    std::string m_wall_het_model{"none"};
    amrex::Real m_monin_obukhov_length{constants::LARGE_NUM};
    amrex::Real m_gamma_m{5.0_rt};
    amrex::Real m_beta_m{16.0_rt};

    //! ``TerrainDrag`` (default) or ``ImmersedTerrain``, read from
    //! ``KLAxell.terrain_model``
    std::string m_terrain_model{"TerrainDrag"};
    //! Wall-model and drag settings shared with ImmersedDragForcing
    std::string m_ib_wall_model{"cell_offset"};
    std::string m_ib_reference_distance{"nominal"};
    amrex::Real m_ib_min_z0{1.0e-4_rt};
    amrex::Real m_ib_drag_coefficient{10.0_rt};
    //! Surface-cell settings shared with ImmersedTerrain
    amrex::Real m_ib_solid_threshold{0.5_rt};
    std::string m_ib_drag_weight{"fraction"};
};

} // namespace kynema_sgf::pde::tke

#endif /* KRANSAXELL_H */
