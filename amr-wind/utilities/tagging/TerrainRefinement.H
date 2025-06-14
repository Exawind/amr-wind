#ifndef TERRAINREFINEMENT_H
#define TERRAINREFINEMENT_H

#include "amr-wind/CFDSim.H"
#include "amr-wind/utilities/tagging/RefinementCriteria.H"
#include "amr-wind/utilities/tagging/Polygon.H"
#include <vector>
#include <sstream>
#include <string>

namespace amr_wind::format_utils {

template <typename T>
inline std::string list_to_string(
    const std::vector<T>& vec,
    const std::string& delimiter = ", ",
    const bool use_parentheses = false)
{
    std::ostringstream oss;
    if (use_parentheses) {
        oss << "(";
    }
    for (size_t i = 0; i < static_cast<size_t>(vec.size()); ++i) {
        if (i > 0) {
            oss << delimiter;
        }
        oss << vec[i];
    }
    if (use_parentheses) {
        oss << ")";
    }
    return oss.str();
}

} // namespace amr_wind::format_utils

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
 * tagging.terrain_lvl1.max_level = 1
 * tagging.terrain_lvl1.grid_buffer_ratio_lo = 0.0 0.5
 * tagging.terrain_lvl1.grid_buffer_ratio_hi = 0.0 0.5
 * tagging.terrain_lvl1.verbose = 0 1
 * tagging.terrain_lvl1.poly_exterior = 0 0  0 10  10 10  10 0  0 0
 * tagging.terrain_lvl1.poly_num_holes = 1
 * tagging.terrain_lvl1.poly_hole_0 = 2 2  2 8  8 8  8 2  2 2
 * tagging.terrain_lvl1.box_lo = 0 0 0
 * tagging.terrain_lvl1.box_hi = 10 10 10
 * \endcode
 *
 * - Only cells within the polygon (if specified), tagging box, and vertical
 *   distance above terrain are refined.
 * - If no polygon is specified, the refinement applies to the whole tagging
 * box.
 * - The grid_buffer_ratio_lo and grid_buffer_ratio_hi parameters control the
 * buffer region below and above the terrain, respectively, as a fraction of the
 * local tile height. Making these ratios larger ensures more cells are tagged
 * for refinement.
 *
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

    //! Terrain file
    std::string m_terrain_file{"terrain.amrwind"};

    //! Ratios of buffer size below the terrain (as a fraction of tile height)
    amrex::Vector<amrex::Real> m_grid_buffer_ratio_lo;

    //! Ratios of buffer size above the terrain (as a fraction of tile height)
    amrex::Vector<amrex::Real> m_grid_buffer_ratio_hi;

    //! Key name for this refinement criteria
    std::string m_key;

    //! Only refine this level if set >= 0
    int m_set_level{-1};

    //! Minimum level to refine (if m_set_level < 0)
    int m_min_level{0};

    //! Maximum level to refine (if m_set_level < 0)
    int m_max_level{32};

    //! Distance above the terrain to refine
    amrex::Real m_vertical_distance{0.0};

    //! The usual tagging bbox */
    amrex::RealBox m_tagging_box;

    //! Optional: The coordinates of the polygon's exterior ring
    amr_wind::polygon_utils::Polygon m_polygon;

    //! Verbose output levels
    amrex::Vector<int> m_verbose_levels;

#ifdef AMREX_USE_GPU
    amrex::Gpu::DeviceVector<amr_wind::polygon_utils::Polygon::Point>
        m_poly_points_dv;
    amrex::Gpu::DeviceVector<int> m_ring_offsets_dv;
#endif

    //! Get active levels as a vector of integers
    AMREX_FORCE_INLINE amrex::Vector<int> get_active_levels() const;

    //! Check if the level should be tagged
    AMREX_FORCE_INLINE bool should_tag_level(int level) const
    {
        if (m_set_level >= 0) {
            return (level == m_set_level);
        }
        return (level >= m_min_level && level <= m_max_level);
    }

    //! Check if the level is verbose
    AMREX_FORCE_INLINE bool is_level_verbose(int level) const
    {
        if (m_verbose_levels.empty()) {
            return false;
        }
        int v = (level < static_cast<int>(m_verbose_levels.size()))
                    ? m_verbose_levels[level]
                    : m_verbose_levels.back();
        return v > 0;
    }
};

} // namespace amr_wind

#endif /* TERRAINREFINEMENT_H */
