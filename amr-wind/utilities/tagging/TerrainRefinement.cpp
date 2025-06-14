#include "amr-wind/CFDSim.H"
#include "AMReX.H"
#include "AMReX_Gpu.H"
#include "AMReX_ParmParse.H"
#include "amr-wind/physics/TerrainDrag.H"
#include "amr-wind/utilities/tagging/TerrainRefinement.H"
#include "amr-wind/utilities/io_utils.H"
#include "amr-wind/utilities/linear_interpolation.H"
#include "amr-wind/utilities/tagging/TerrainElevation.H"

namespace amr_wind {

TerrainRefinement::TerrainRefinement(const CFDSim& sim)
    : m_sim(sim), m_tagging_box(m_sim.repo().mesh().Geom(0).ProbDomain())
{}

void TerrainRefinement::initialize(const std::string& key)
{
    m_key = key.substr(8); // Store the key name
    amrex::ParmParse pp(key);
    pp.queryarr("verbose", m_verbose_levels);
    pp.queryarr("grid_buffer_ratio_lo", m_grid_buffer_ratio_lo);
    pp.queryarr("grid_buffer_ratio_hi", m_grid_buffer_ratio_hi);

    //! Reading the Terrain Coordinates from  file
    amrex::ParmParse pp_terrain("TerrainDrag");
    pp_terrain.query("terrain_file", m_terrain_file);
    amr_wind::TerrainElevation::ensure_loaded(m_terrain_file);

    // Outer radial extent of the cylinder, always read in from input file
    pp.get("vertical_distance", m_vertical_distance);

    // Get the levels
    pp.query("level", m_set_level);
    if (m_set_level >= 0) {
        m_min_level = m_set_level;
        m_max_level = m_set_level;
    } else {
        pp.query("min_level", m_min_level);
        pp.query("max_level", m_max_level);
    }

    // Read polygon from input file
    m_polygon.read_from_parmparse(key);

    amrex::Vector<amrex::Real> box_lo(AMREX_SPACEDIM, 0);
    amrex::Vector<amrex::Real> box_hi(AMREX_SPACEDIM, 0);
    if (pp.queryarr("box_lo", box_lo, 0, static_cast<int>(box_lo.size())) ==
        1) {
        m_tagging_box.setLo(box_lo);
    }
    if (pp.queryarr("box_hi", box_hi, 0, static_cast<int>(box_hi.size())) ==
        1) {
        m_tagging_box.setHi(box_hi);
    }

#ifdef AMREX_USE_GPU
    m_poly_points_dv.resize(m_polygon.points().size());
    amrex::Gpu::copyAsync(
        amrex::Gpu::hostToDevice, m_polygon.points().begin(),
        m_polygon.points().end(), m_poly_points_dv.begin());

    m_ring_offsets_dv.resize(m_polygon.ring_offsets().size());
    amrex::Gpu::copyAsync(
        amrex::Gpu::hostToDevice, m_polygon.ring_offsets().begin(),
        m_polygon.ring_offsets().end(), m_ring_offsets_dv.begin());
#endif

    auto active_levels =
        amr_wind::format_utils::list_to_string(get_active_levels(), ", ", true);

    amrex::Print() << "Created terrain refinement " << m_key << " for levels "
                   << active_levels << " and vertical distance "
                   << m_vertical_distance << " (grid_buffer_ratio_lo: "
                   << amr_wind::format_utils::list_to_string(
                          m_grid_buffer_ratio_lo, ", ", true)
                   << ", grid_buffer_ratio_hi: "
                   << amr_wind::format_utils::list_to_string(
                          m_grid_buffer_ratio_hi, ", ", true)
                   << ")" << std::endl;
}

void TerrainRefinement::operator()(
    int level, amrex::TagBoxArray& tags, amrex::Real /*time*/, int /*ngrow*/)
{
    const bool verbose = is_level_verbose(level);

    if (!should_tag_level(level)) {
        return;
    }

    // Get the buffer ratios for this level
    auto buffer_lo = m_grid_buffer_ratio_lo.empty()
                         ? 0.0
                         : (level < m_grid_buffer_ratio_lo.size()
                                ? m_grid_buffer_ratio_lo[level]
                                : m_grid_buffer_ratio_lo.back());
    auto buffer_hi = m_grid_buffer_ratio_hi.empty()
                         ? 0.0
                         : (level < m_grid_buffer_ratio_hi.size()
                                ? m_grid_buffer_ratio_hi[level]
                                : m_grid_buffer_ratio_hi.back());

    // Geometry
    const auto& geom = m_sim.repo().mesh().Geom(level);
    const auto& prob_lo = geom.ProbLoArray();
    const auto& dx = geom.CellSizeArray();

    // Get tagging box and vertical distance for device
    const auto tagging_box = m_tagging_box;

    // Round up to the nearest positive multiple of dx[2] (if not already a
    // multiple)
    auto vertical_distance = m_vertical_distance;
    if (vertical_distance <= 0.0) {
        vertical_distance = dx[2];
    } else {
        vertical_distance = std::ceil(vertical_distance / dx[2]) * dx[2];
    }

    // Tag arrays
    const auto& tag_arrs = tags.arrays();

    const auto& mfab = (*m_sim.repo().fields()[0])(level);

    // Prepare polygon data for GPU
    const bool polygon_is_empty = m_polygon.is_empty();
    auto n_rings = m_polygon.num_rings();
    auto n_points = m_polygon.points().size();

#ifdef AMREX_USE_GPU
    auto const* p_poly_points = m_poly_points_dv.data();
    auto const* p_ring_offsets = m_ring_offsets_dv.data();
#else
    auto const* p_poly_points = m_polygon.points().data();
    auto const* p_ring_offsets = m_polygon.ring_offsets().data();
#endif

    // Get polygon bounding box (lo, hi)
    amr_wind::polygon_utils::Polygon::Point bbox_lo = {0.0, 0.0};
    amr_wind::polygon_utils::Polygon::Point bbox_hi = {0.0, 0.0};
    if (!polygon_is_empty) {
        if (!m_polygon.bbox_valid()) {
            m_polygon.compute_bounding_box();
        }
        m_polygon.get_bounding_box(bbox_lo, bbox_hi);
    }

    auto terrain_ptr = amr_wind::TerrainElevation::get_instance();
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        terrain_ptr && terrain_ptr->is_loaded(),
        "TerrainElevation must be loaded before TerrainRefinement is used.");

    const auto& xterrain = terrain_ptr->x();
    const auto& yterrain = terrain_ptr->y();
    const auto& zterrain = terrain_ptr->z();

    const auto xterrain_size = xterrain.size();
    const auto yterrain_size = yterrain.size();
    const auto zterrain_size = zterrain.size();
    amrex::Gpu::DeviceVector<amrex::Real> d_xterrain(xterrain_size);
    amrex::Gpu::DeviceVector<amrex::Real> d_yterrain(yterrain_size);
    amrex::Gpu::DeviceVector<amrex::Real> d_zterrain(zterrain_size);
    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, xterrain.begin(), xterrain.end(),
        d_xterrain.begin());
    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, yterrain.begin(), yterrain.end(),
        d_yterrain.begin());
    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, zterrain.begin(), zterrain.end(),
        d_zterrain.begin());
    const auto* xterrain_ptr = d_xterrain.data();
    const auto* yterrain_ptr = d_yterrain.data();
    const auto* zterrain_ptr = d_zterrain.data();

    for (amrex::MFIter mfi(mfab, amrex::TilingIfNotGPU()); mfi.isValid();
         ++mfi) {
        const auto& bx = mfi.tilebox();
        int nbx = mfi.LocalIndex();

        // Compute buffer height based on the vertical distance and buffer ratio
        // This buffer is used to ensure that we tag cells that are slightly
        // above the terrain, allowing for a smoother transition in terrain
        // refinement.
        // The buffer is a fraction of the box height in the z-direction
        // and is defined as a ratio of the vertical distance.
        // This ensures that the buffer is proportional to the vertical distance
        // and adapts to the size of the box in the z-direction.
        const auto box_height = bx.length(2) * dx[2];
        auto buffer_z_lo = buffer_lo * box_height;
        buffer_z_lo = std::max(buffer_z_lo, 0.5 * dx[2]);
        auto buffer_z_hi = buffer_hi * box_height;

        // Compute tilebox bounds in physical space
        auto tile_lo_x = prob_lo[0] + (bx.smallEnd(0) + 0.0) * dx[0];
        auto tile_hi_x = prob_lo[0] + (bx.bigEnd(0) + 1.0) * dx[0];
        auto tile_lo_y = prob_lo[1] + (bx.smallEnd(1) + 0.0) * dx[1];
        auto tile_hi_y = prob_lo[1] + (bx.bigEnd(1) + 1.0) * dx[1];

        // Only check polygon if not empty
        if (!polygon_is_empty) {
            // Check for overlap with polygon bounding box
            if (tile_hi_x < bbox_lo[0] || tile_lo_x > bbox_hi[0] ||
                tile_hi_y < bbox_lo[1] || tile_lo_y > bbox_hi[1]) {
                if (verbose) {
                    amrex::Print()
                        << m_key << " - Level " << level
                        << " (vertical_distance=" << vertical_distance
                        << "): Skipping tile: [" << tile_lo_x << ","
                        << tile_hi_x << "] x [" << tile_lo_y << "," << tile_hi_y
                        << "] "
                        << "does not overlap polygon bbox [" << bbox_lo[0]
                        << "," << bbox_hi[0] << "] x [" << bbox_lo[1] << ","
                        << bbox_hi[1] << "]\n";
                }
                continue; // Skip this tile, no overlap
            }
        }

        amrex::ParallelFor(
            bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                // 1. Early exit if already tagged
                if (tag_arrs[nbx](i, j, k) == amrex::TagBox::SET) {
                    return;
                }

                const amrex::RealVect coord = {AMREX_D_DECL(
                    prob_lo[0] + (i + 0.5) * dx[0],
                    prob_lo[1] + (j + 0.5) * dx[1],
                    prob_lo[2] + (k + 0.5) * dx[2])};

                const amrex::Real terrainHt = interp::bilinear(
                    xterrain_ptr, xterrain_ptr + xterrain_size, yterrain_ptr,
                    yterrain_ptr + yterrain_size, zterrain_ptr, coord[0],
                    coord[1]);

                const auto cellHt = coord[2] - terrainHt;

                const amr_wind::polygon_utils::Polygon::Point center_point{
                    coord[0], coord[1]};

                // 2. Check vertical distance
                if ((cellHt < -buffer_z_lo) ||
                    (cellHt > std::max(vertical_distance, buffer_z_hi))) {
                    return;
                }

                // 3. Check tagging box
                if (!tagging_box.contains(coord)) {
                    return;
                }

                // 4. Check polygon bounding box
                if (!polygon_is_empty) {
                    if (!amr_wind::polygon_utils::Polygon::
                            poly_bounding_box_contains(
                                center_point, bbox_lo, bbox_hi)) {
                        return;
                    }
                }

                // 5. Point-in-polygon
                bool in_poly = false;
                if (!polygon_is_empty) {
                    in_poly =
                        amr_wind::polygon_utils::Polygon::is_point_in_polygon(
                            p_poly_points, p_ring_offsets, n_rings,
                            static_cast<int>(n_points), center_point);

                } else {
                    in_poly = true;
                }

                if (in_poly) {
                    tag_arrs[nbx](i, j, k) = amrex::TagBox::SET;
                }
            });
    }

    amrex::Gpu::streamSynchronize();
}

amrex::Vector<int> TerrainRefinement::get_active_levels() const
{
    const int max_level = m_sim.repo().mesh().maxLevel();
    amrex::Vector<int> levels;
    int min_lev = std::max(0, m_set_level >= 0 ? m_set_level : m_min_level);
    int max_lev =
        std::min(max_level, m_set_level >= 0 ? m_set_level : m_max_level);

    for (int lev = min_lev; lev <= max_lev; ++lev) {
        levels.push_back(lev);
    }
    return levels;
}

} // namespace amr_wind