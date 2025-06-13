#ifndef TERRAINREFINEMENT_H
#define TERRAINREFINEMENT_H

#include "amr-wind/CFDSim.H"
#include "amr-wind/utilities/tagging/RefinementCriteria.H"
#include "amr-wind/utilities/tagging/Polygon.H"

namespace amr_wind {
class Field;

/**
 * \class TerrainRefinement
 * \ingroup amr_utils
 * \brief AMR refinement using a given distance to terrain and optional polygon
 * region.
 *
 * Example usage in input file:
 * \code
 * tagging.labels = terrain_lvl1
 * tagging.terrain_lvl1.type = TerrainRefinement
 * tagging.terrain_lvl1.vertical_distance = 100
 * tagging.terrain_lvl1.level = 0
 * tagging.terrain_lvl1.poly_exterior = 0 0  0 10  10 10  10 0  0 0
 * tagging.terrain_lvl1.poly_num_holes = 1
 * tagging.terrain_lvl1.poly_hole_0 = 2 2  2 8  8 8  8 2  2 2
 * tagging.terrain_lvl1.box_lo = 0 0 0
 * tagging.terrain_lvl1.box_hi = 10 10 10
 * \endcode
 *
 * - Only cells within the polygon (if specified), tagging box, and vertical
 * distance above terrain are refined.
 * - If no polygon is specified, the refinement applies to the whole tagging
 * box.
 */

class TerrainRefinement : public RefinementCriteria::Register<TerrainRefinement>
{
public:
    static std::string identifier() { return "TerrainRefinement"; }

    explicit TerrainRefinement(const CFDSim& sim);

    ~TerrainRefinement() override = default;

    /**
     * \brief Read input file and initialize boxarray used to refine each level
     * \param key ParmParse prefix for this refinement
     */
    void initialize(const std::string& key) override;

    /**
     * \brief Tag cells for refinement
     * \param level AMR level
     * \param tags TagBoxArray to set
     * \param time Current simulation time
     * \param ngrow Number of ghost cells
     */
    void operator()(
        const int level,
        amrex::TagBoxArray& tags,
        const amrex::Real time,
        const int ngrow) override;

private:
    const CFDSim& m_sim;
    /** \brief Key name for this refinement criteria */
    std::string m_key;
    /** \brief Pointer to the terrain_height field */
    Field* m_terrain_height{nullptr};
    /** \brief Pointer to the terrain_blank field so that we try not to refine
     * under the terrain */
    IntField* m_terrain_blank{nullptr};
    /** \brief Grid refinement level */
    int m_max_lev{0};
    /** \brief Distance above the terrain to refine */
    amrex::Real m_vertical_distance{0.0};
    /** \brief The usual tagging bbox */
    amrex::RealBox m_tagging_box;
    /** \brief Optional: The coordinates of the polygon's exterior ring */
    amr_wind::polygon_utils::Polygon m_polygon;
    /** \brief Verbose output level */
    int m_verbose{-1};

#ifdef AMREX_USE_GPU
    amrex::Gpu::DeviceVector<amr_wind::polygon_utils::Polygon::Point>
        m_poly_points_dv;
    amrex::Gpu::DeviceVector<int> m_ring_offsets_dv;
#endif
};

} // namespace amr_wind

#endif /* TERRAINREFINEMENT_H */
