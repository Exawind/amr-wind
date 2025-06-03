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
    if (m_max_lev <= 0) {
        amrex::Abort("TerrainRefinement: level must be > 0");
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

    amrex::Print() << "Created terrain refinement with level " << m_max_lev
                   << " and vertical distance " << m_vertical_distance
                   << std::endl;
}

void TerrainRefinement::operator()(
    int level, amrex::TagBoxArray& tags, amrex::Real time, int ngrow)
{
    (void)time;
    (void)ngrow;
    const bool do_tag = level < m_max_lev;
    if (!do_tag) {
        return;
    }
    // Geometry
    const auto& geom = m_sim.repo().mesh().Geom(level);
    const auto& prob_lo = geom.ProbLoArray();
    const auto& dx = geom.CellSizeArray();

    // Get tagging box and vertical distance for device
    const auto tagging_box = m_tagging_box;
    auto vertical_distance = m_vertical_distance;

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

    amrex::ParallelFor(
        mfab, [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
            const amrex::Real z = prob_lo[2] + (k + 0.5) * dx[2];
            const amrex::Real terrainHt = mterrain_h_arrs[nbx](i, j, k);
            const amrex::Real cellHt = z - terrainHt;

            const amrex::RealVect coord = {AMREX_D_DECL(
                prob_lo[0] + (i + 0.5) * dx[0], prob_lo[1] + (j + 0.5) * dx[1],
                prob_lo[2] + (k + 0.5) * dx[2])};

            const amr_wind::polygon_utils::Polygon::Point testPt{
                coord[0], coord[1]};

            // 1. Check tagging box first
            if (!tagging_box.contains(coord)) {
                return;
            }

            // 2. Check vertical distance
            if ((cellHt < -0.5 * dx[2]) || (cellHt > vertical_distance)) {
                return;
            }

            // 3. Check terrain blanking
            if (mterrain_b_arrs[nbx](i, j, k) >= 1) {
                return;
            }

            // 4. Polygon point-in-polygon test using the new layout
            bool in_poly = false;
            if (!polygon_is_empty) {
                // Outer ring
                int start = static_cast<int>(p_ring_offsets[0]);
                int end = (n_rings > 1) ? static_cast<int>(p_ring_offsets[1])
                                        : static_cast<int>(n_points);
                int n = end - start;
                if (amr_wind::polygon_utils::Polygon::is_point_in_ring(
                        p_poly_points + start, n, testPt)) {
                    in_poly = true;
                    // Check holes
                    for (int ring_i = 1; ring_i < n_rings; ++ring_i) {
                        int h_start = static_cast<int>(p_ring_offsets[ring_i]);
                        int h_end =
                            (ring_i + 1 < n_rings)
                                ? static_cast<int>(p_ring_offsets[ring_i + 1])
                                : static_cast<int>(n_points);
                        int h_n = h_end - h_start;
                        if (amr_wind::polygon_utils::Polygon::is_point_in_ring(
                                p_poly_points + h_start, h_n, testPt)) {
                            in_poly = false;
                            break;
                        }
                    }
                }
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