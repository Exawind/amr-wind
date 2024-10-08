#ifndef TerrainDrag_H
#define TerrainDrag_H

#include "amr-wind/core/Physics.H"
#include "amr-wind/core/Field.H"
#include "amr-wind/CFDSim.H"

namespace amr_wind::terraindrag {

namespace {} // namespace

/** Terraindrag Flow physics
 *  \ingroup physics
 */

class TerrainDrag : public Physics::Register<TerrainDrag>
{
public:
    static std::string identifier() { return "TerrainDrag"; }

    explicit TerrainDrag(CFDSim& sim);

    ~TerrainDrag() override = default;

    void
    initialize_fields(int /*level*/, const amrex::Geometry& /*geom*/) override
    {}

    void pre_init_actions() override;

    void post_init_actions() override;

    void post_regrid_actions() override {}

    void pre_advance_work() override {}

    void post_advance_work() override {}

    int return_blank_value(int i, int j, int k);

private:
    CFDSim& m_sim;
    const FieldRepo& m_repo;
    const amrex::AmrCore& m_mesh;
    Field& m_velocity;
    //! Blanking Field for Terrain or Buildings
    IntField& m_terrain_blank;
    //! Terrain Drag Force Term
    IntField& m_terrain_drag;
    //! Reading the Terrain Coordinates from  file
    amrex::Vector<amrex::Real> m_xterrain;
    amrex::Vector<amrex::Real> m_yterrain;
    amrex::Vector<amrex::Real> m_zterrain;
    //! Roughness Field
    Field& m_terrainz0;
    Field& m_terrain_height;
    //! Reading the Roughness Coordinates from file - Not Fully there yet
    //! Need updates to ABLWallFunction in future
    amrex::Vector<amrex::Real> m_xrough;
    amrex::Vector<amrex::Real> m_yrough;
    amrex::Vector<amrex::Real> m_z0rough;
};
} // namespace amr_wind::terraindrag

#endif
