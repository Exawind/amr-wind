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

    void initialize_fields(int level, const amrex::Geometry& geom) override;

    void post_init_actions() override {}

    void post_regrid_actions() override {}

    void pre_advance_work() override {}

    void post_advance_work() override {}

    //    amrex::Real computeTerrainHeight(const amrex::Real x1,const
    //    amrex::Real x2) const;
private:
    CFDSim& m_sim;
    const FieldRepo& m_repo;
    const amrex::AmrCore& m_mesh;
    Field& m_velocity;
    // Blanking Field for Terrain or Buildings
    Field& m_terrainBlank;
    // Terrain Height field - Not in use now
    Field& m_terrainDrag;
    // Reading the Terrain Coordinates from  file
    amrex::Vector<amrex::Real> m_xterrain, m_yterrain, m_zterrain;
};
} // namespace amr_wind::terraindrag

#endif
