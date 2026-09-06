#ifndef IMMERSEDTERRAIN_H
#define IMMERSEDTERRAIN_H

#include "src/core/Physics.H"
#include "src/core/Field.H"
#include "src/CFDSim.H"

namespace kynema_sgf::immersedterrain {

/** Immersed terrain geometry for the drag-based immersed boundary method
 *  \ingroup physics
 *
 *  Successor to TerrainDrag with a partial (0..1) terrain volume fraction
 *  instead of a binary blanking, so that slopes are not represented as a
 *  staircase. The class only builds geometry fields; the momentum and
 *  temperature forcing consume them.
 *
 *  Fields declared (all with one ghost cell):
 *
 *  - ``terrain_fraction`` (Real, 1 component): fraction of the cell volume
 *    occupied by terrain. 0 in fluid, 1 inside terrain, in between for cells
 *    cut by the terrain surface. The immersed boundary drag is scaled by it.
 *
 *  - ``terrain_mask`` (Int, 1 component): cell classification
 *      - 0 ``mask_fluid``   : fluid cell away from the terrain
 *      - 1 ``mask_solid``   : cell fully inside the terrain
 *      - 2 ``mask_surface`` : cell where the wall model (MOST) is applied:
 *        a partially filled cell, or a fluid cell sharing a face with a
 *        mostly solid cell on any of its six sides. The six-way search is
 *        needed to catch the vertical faces of steep terrain and buildings
 *        that a purely vertical search misses.
 *    Use with ``FieldRefinement``: ``field_error = 1.5`` tags only the
 *    surface band, ``field_error = 0.5`` tags surface and solid cells.
 *
 *  - ``terrain_surface`` (Real, 3 components): terrain height h(x,y) at the
 *    cell center, and the slopes dh/dx and dh/dy. The slopes are what a
 *    second-order (planar) reconstruction of the surface inside a cell needs;
 *    the current fraction uses only the height and is first order.
 *
 *  - ``terrain_roughness`` (Real, 1 component): aerodynamic roughness z0,
 *    uniform or interpolated from a roughness file.
 *
 *  - ``terrain_drag_rate`` (Real, 1 component), only with
 *    ``implicit_projection = true``: beta * C_d / dz, the immersed drag rate
 *    of each cell. Its presence makes the nodal and MAC projections treat
 *    the terrain as a fluid of density rho (1 + rate dt), so that the drag is
 *    applied implicitly with the pressure and the pressure gradient produces
 *    no velocity inside the terrain. ImmersedDragForcing then skips its
 *    explicit drag term. Without this, the projection re-injects dt grad(p)/rho
 *    inside the body every step and the interior residual is first order in dt.
 *
 *  - ``terrain_diffusion_xf/yf/zf`` (Real, face centered), only with
 *    ``interface_diffusion`` other than ``none``: factor applied to the face
 *    diffusion coefficients of every equation solved with the shared diffusion
 *    operator. ``block`` uses min(1 - beta_L, 1 - beta_R), which removes the
 *    molecular and SGS flux across faces touching the terrain so that the wall
 *    model alone carries the wall stress (turbulent pathway). ``no_slip``
 *    scales the coefficient of fluid/solid faces by dx_f / d1 so that the
 *    discrete flux mu (u - 0)/dx_f equals the flux mu u / d1 to a no-slip wall
 *    at the true wall position (laminar pathway); d1 = z - h on the bottom
 *    face, dx_f / 2 elsewhere.
 */
class ImmersedTerrain : public Physics::Register<ImmersedTerrain>
{
public:
    static std::string identifier() { return "ImmersedTerrain"; }

    //! Values stored in the ``terrain_mask`` int field
    static constexpr int mask_fluid = 0;
    static constexpr int mask_solid = 1;
    static constexpr int mask_surface = 2;

    //! Component indices of the ``terrain_surface`` field
    static constexpr int surf_height = 0;
    static constexpr int surf_slope_x = 1;
    static constexpr int surf_slope_y = 2;

    explicit ImmersedTerrain(CFDSim& sim);

    ~ImmersedTerrain() override = default;

    void initialize_fields(int level, const amrex::Geometry& geom) override;

    void pre_init_actions() override {}

    void post_init_actions() override {}

    void post_regrid_actions() override;

    void pre_advance_work() override {}

    void post_advance_work() override {}

private:
    CFDSim& m_sim;
    const FieldRepo& m_repo;
    const amrex::AmrCore& m_mesh;

    //! Terrain volume fraction per cell (0..1)
    Field& m_terrain_fraction;
    //! Cell classification: fluid / solid / surface
    IntField& m_terrain_mask;
    //! Terrain height and slopes at the cell center
    Field& m_terrain_surface;
    //! Aerodynamic roughness
    Field& m_terrain_roughness;

    //! Terrain file (x, y, h on a flat grid)
    std::string m_terrain_file{"terrain.amrwind"};
    //! Optional roughness file (x, y, z0 on a flat grid)
    std::string m_roughness_file{"terrain.roughness"};
    //! Roughness used when no roughness file is present
    amrex::Real m_uniform_z0{0.1_rt};

    //! ``volume_fraction`` (column fraction of the cell below the terrain) or
    //! ``distance_function`` (tanh of signed distance / smoothing length)
    std::string m_blanking_method{"volume_fraction"};
    //! Smoothing length for ``distance_function`` in units of dz
    amrex::Real m_smoothing_length{1.0_rt};
    //! Neighbor is treated as solid for the surface search above this
    amrex::Real m_solid_threshold{0.5_rt};

    //! Apply the immersed drag implicitly through the projections
    bool m_implicit_projection{false};
    //! Drag coefficient C_d shared with ImmersedDragForcing
    amrex::Real m_drag_coefficient{10.0_rt};
    //! beta * C_d / dz, declared only when m_implicit_projection is true
    Field* m_terrain_drag_rate{nullptr};

    //! ``none``, ``block`` or ``no_slip``
    std::string m_interface_diffusion{"none"};
    //! Face factors for the diffusion coefficients, declared only when
    //! m_interface_diffusion is not ``none``
    amrex::Array<Field*, AMREX_SPACEDIM> m_diffusion_factor{
        {nullptr, nullptr, nullptr}};
};

} // namespace kynema_sgf::immersedterrain

#endif
