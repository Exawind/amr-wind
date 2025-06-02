#ifndef FIELDREFINEMENT_H
#define FIELDREFINEMENT_H

#include "amr-wind/CFDSim.H"
#include "amr-wind/utilities/tagging/poly_ops.H"
#include "amr-wind/utilities/tagging/RefinementCriteria.H"

namespace amr_wind {
// class CFDSim;
class Field;
// class IntField;

/** AMR refinement using a given distance to terrain
 *  \ingroup amr_utils
 *
 *  ```
 *  tagging.labels = terrain_lvl1 terrain_lvl2
 *  tagging.terrain_lvl1.type = TerrainRefinement
 *  tagging.terrain_lvl1.vertical_distance = 100
 *  tagging.terrain_lvl1.level = 1
 *  tagging.terrain_lvl1.poly_outer = ...
 *  tagging.terrain_lvl1.poly_inners = ...
 *  tagging.terrain_lvl1.box_lo = ...
 *  tagging.terrain_lvl1.box_hi = ...
 *  tagging.terrain_lvl2.type = TerrainRefinement
 *  tagging.terrain_lvl2.vertical_distance = 40
 *  tagging.terrain_lvl2.level = 2
 *  tagging.terrain_lvl2.poly_outer = ...
 *
 * poly_outer should be define as: number_of_points pt0_x pt0_y ... ptN_x ptN_y
 * poly_outer poly_inners be define as: number_of_inner_rings
 * number_of_points_in_ring0 r0_pt0_x r0_pt0_x ... r0_ptN_x r0_ptN_x
 * number_of_points_in_ringN rN_pt0_x rN_pt0_x ... rN_ptN_x rN_ptN_x
 *  ```
 */

class TerrainRefinement : public RefinementCriteria::Register<TerrainRefinement>
{
public:
    static std::string identifier() { return "TerrainRefinement"; }

    explicit TerrainRefinement(const CFDSim& sim);

    ~TerrainRefinement() override = default;

    //! Read input file and initialize boxarray used to refine each level
    void initialize(const std::string& key) override;

    void operator()(
        const int level,
        amrex::TagBoxArray& tags,
        const amrex::Real time,
        const int ngrow) override;

private:
    const CFDSim& m_sim;

    // Pointer to the terrain_height field.
    Field* m_terrain_height{nullptr};
    // Pointer to the terrain_blank field so that we try not to refine under the
    // terrain
    IntField* m_terrain_blank{nullptr};
    // Target grid refinement level
    int m_max_lev{0};
    // Distance above the terrain to refine
    amrex::Real m_vertical_distance{0.0};
    // The usual tagging bbox
    amrex::RealBox m_tagging_box;

    // Optional: The coordinates of the polygon's exterior ring
    amrex::Vector<amr_wind::polygon_utils::Point> m_poly_outer;
    // Optional: The coordinates of the polygon's interior ring(s)
    amrex::Vector<amrex::Vector<amr_wind::polygon_utils::Point>> m_poly_rings;
};

} // namespace amr_wind

#endif /* TERRAINREFINEMENT_H */
