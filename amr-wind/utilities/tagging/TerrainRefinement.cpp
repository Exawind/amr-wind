#include "amr-wind/utilities/tagging/TerrainRefinement.H"
#include "amr-wind/CFDSim.H"

#include "AMReX.H"
#include "AMReX_ParmParse.H"

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
        amrex::Abort("TerrainRefinement: level should be strictly above 0");
    }

    amrex::Vector<amrex::Real> poly_outer;
    amrex::Vector<amrex::Real> poly_inners;
    pp.queryarr("poly_outer", poly_outer);

    if (poly_outer.size() > 0) {            
        pp.getarr("poly_outer", poly_outer);
        const int n = poly_outer.size();
        const int n_points = static_cast<int>(poly_outer[0]);
        // amrex::Print() << "Number of points:" << n_points << std::endl;
        const int n_expected = n_points * 2 + 1;
        if (n_expected != n) {
            amrex::Abort(
                "Expected a list of " + std::to_string(n_expected) +
                " numbers, found " + std::to_string(n - 1) + "!");
        }
        m_poly_outer.resize(n_points);
        for (int i = 0; i < n_points; ++i) {
            const auto pt_x = poly_outer[1 + 2 * i];
            const auto pt_y = poly_outer[1 + 2 * i + 1];
            m_poly_outer[i] = amr_wind::polygon_utils::Point({pt_x, pt_y});
            // amrex::Print() << pt_x << " " << pt_y << std::endl;
        }
        pp.queryarr("poly_inners", poly_inners);
        if (poly_inners.size() > 0) {            
            const int n_rings = static_cast<int>(poly_inners[0]);
            m_poly_rings.resize(n_rings);
            // amrex::Print() << "Number of inner rings:" << n_rings <<
            // std::endl;
            int offset = 1;
            for (int ring_i = 0; ring_i < n_rings; ++ring_i) {
                const int n_pts = static_cast<int>(poly_inners[offset]);
                m_poly_rings[ring_i].resize(n_pts);
                // amrex::Print() << "Number of inner points:" << n_pts <<
                // std::endl;
                offset += 1;
                for (int pt_i = 0; pt_i < n_pts; ++pt_i) {
                    const auto pt_x = poly_inners[offset + 2 * pt_i];
                    const auto pt_y = poly_inners[offset + 2 * pt_i + 1];
                    m_poly_rings[ring_i][pt_i] =
                        amr_wind::polygon_utils::Point({pt_x, pt_y});
                    // amrex::Print() << m_poly_rings[ring_i][pt_i].x <<
                    // std::endl;
                }
                offset += 2 * n_pts;
            }
            if ((n_rings > 0) && (offset != poly_inners.size())) {
                amrex::Abort(
                    "Expected a list of " + std::to_string(offset) +
                    " numbers, found " + std::to_string(poly_inners.size()) +
                    "!");
            }
        }
    }

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
    int level, amrex::TagBoxArray& tags, amrex::Real time, int /*ngrow*/)
{
    const bool do_tag = level < m_max_lev;
    if (!do_tag) {
        return;
    }
    const auto& repo = m_sim.repo();
    const auto& geom = repo.mesh().Geom(level);
    const auto& prob_lo = geom.ProbLoArray();
    const auto& dx = geom.CellSizeArray();
    const auto tagging_box = m_tagging_box;

    const auto& tag_arrs = tags.arrays();

    // const auto& farrs = mfab.const_arrays();
    const auto& mfab = (*m_terrain_height)(level);
    const auto& mterrain_h_arrs = mfab.const_arrays();
    const auto& mterrain_b_arrs = (*m_terrain_blank)(level).const_arrays();

    amrex::ParallelFor(
        mfab, [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
            const amrex::Real z = prob_lo[2] + (k + 0.5) * dx[2];
            const amrex::Real terrainHt = mterrain_h_arrs[nbx](i, j, k);
            const amrex::Real cellHt = z - terrainHt;

            const amrex::RealVect coord = {AMREX_D_DECL(
                prob_lo[0] + (i + 0.5) * dx[0], prob_lo[1] + (j + 0.5) * dx[1],
                prob_lo[2] + (k + 0.5) * dx[2])};

            const auto testPt =
                amr_wind::polygon_utils::Point({coord[0], coord[1]});
            bool in_poly = m_poly_outer.size()
                               ? amr_wind::polygon_utils::is_point_in_polygon(
                                     m_poly_outer, testPt)
                               : true;
            if (in_poly) {
                for (int ring_i = 0; ring_i < m_poly_rings.size(); ++ring_i) {
                    const bool in_ring =
                        amr_wind::polygon_utils::is_point_in_polygon(
                            m_poly_rings[ring_i], testPt);
                    if (in_ring) {
                        in_poly = false;
                        break;
                    }
                }
            }

            if (((cellHt >= -0.5) * dx[2] && (cellHt <= m_vertical_distance)) &&
                (mterrain_b_arrs[nbx](i, j, k) < 1) && in_poly &&
                (tagging_box.contains(coord))) {
                tag_arrs[nbx](i, j, k) = amrex::TagBox::SET;
            }
        });
}

} // namespace amr_wind
