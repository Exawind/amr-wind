#include "amr-wind/CFDSim.H"
#include "AMReX.H"
#include "AMReX_Gpu.H"
#include "AMReX_ParmParse.H"
#include "amr-wind/utilities/tagging/TerrainRefinement.H"

namespace amr_wind {

TerrainRefinement::TerrainRefinement(const CFDSim& sim)
    : m_sim(sim), m_tagging_box(m_sim.repo().mesh().Geom(0).ProbDomain())
{}

void TerrainRefinement::initialize(const std::string& key)
{
    amrex::ParmParse pp(key);

    const auto& repo = m_sim.repo();

    const bool is_terrain = repo.field_exists("terrain_height");
    if (!is_terrain) {
        amrex::Abort("Need terrain blanking variable to use this refinement");
    }
    m_terrain_height = &(m_sim.repo().get_field("terrain_height"));
    m_terrain_blank = &(m_sim.repo().get_int_field("terrain_blank"));

    // Outer radial extent of the cylinder, always read in from input file
    pp.get("vertical_distance", m_vertical_distance);
    pp.get("level", m_max_lev);
    if (m_max_lev < 0) {
        amrex::Abort("TerrainRefinement: level must be >= 0");
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

    amrex::Print() << "Created terrain refinement for level " << m_max_lev
                   << " and vertical distance " << m_vertical_distance
                   << std::endl;
}

void TerrainRefinement::operator()(
    int level, amrex::TagBoxArray& tags, amrex::Real time, int ngrow)
{
    (void)time;
    (void)ngrow;
    const bool do_tag = level <= m_max_lev;
    if (!do_tag) {
        return;
    }
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

    // Get terrain arrays
    const auto& mfab = (*m_terrain_height)(level);
    const auto& mterrain_h_arrs = mfab.const_arrays();
    const auto& mterrain_b_arrs = (*m_terrain_blank)(level).const_arrays();

    // Prepare polygon data for GPU
    const bool polygon_is_empty = m_polygon.is_empty();
    auto n_rings = m_polygon.num_rings();
    auto n_points = m_polygon.points().size();

#ifdef AMREX_USE_GPU
    amrex::Gpu::DeviceVector<amr_wind::polygon_utils::Polygon::Point>
        poly_points_dv(n_points);
    amrex::Gpu::copyAsync(
        amrex::Gpu::hostToDevice, m_polygon.points().begin(),
        m_polygon.points().end(), poly_points_dv.begin());
    auto const* p_poly_points = poly_points_dv.data();

    amrex::Gpu::DeviceVector<int> ring_offsets_dv(
        m_polygon.ring_offsets().size());
    amrex::Gpu::copyAsync(
        amrex::Gpu::hostToDevice, m_polygon.ring_offsets().begin(),
        m_polygon.ring_offsets().end(), ring_offsets_dv.begin());
    auto const* p_ring_offsets = ring_offsets_dv.data();
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

    amrex::ParallelFor(
        mfab, [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
            // 1. Early exit if already tagged
            if (tag_arrs[nbx](i, j, k) == amrex::TagBox::SET) {
                return;
            }

            const amrex::RealVect coord = {AMREX_D_DECL(
                prob_lo[0] + (i + 0.5) * dx[0], prob_lo[1] + (j + 0.5) * dx[1],
                prob_lo[2] + (k + 0.5) * dx[2])};

            const auto terrainHt = mterrain_h_arrs[nbx](i, j, k);
            const auto cellHt = coord[2] - terrainHt;

            const amr_wind::polygon_utils::Polygon::Point center_point{
                coord[0], coord[1]};

            // 1. Check terrain blanking
            if (mterrain_b_arrs[nbx](i, j, k) >= 1) {
                return;
            }

            // 2. Check vertical distance
            if ((cellHt < -0.5 * dx[2]) || (cellHt > vertical_distance)) {
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
                in_poly = amr_wind::polygon_utils::Polygon::is_point_in_polygon(
                    p_poly_points, p_ring_offsets, n_rings,
                    static_cast<int>(n_points), center_point);

            } else {
                in_poly = true;
            }

            if (in_poly) {
                tag_arrs[nbx](i, j, k) = amrex::TagBox::SET;
            }
        });

    amrex::Gpu::streamSynchronize();
}

} // namespace amr_wind