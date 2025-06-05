#ifndef ForestDrag_H
#define ForestDrag_H

#include "amr-wind/core/Physics.H"
#include "amr-wind/core/Field.H"
#include "amr-wind/CFDSim.H"
#include "amr-wind/utilities/index_operations.H"
#include "amr-wind/utilities/integrals.H"
#include "amr-wind/utilities/constants.H"

namespace amr_wind::forestdrag {

struct Forest
{
    int m_id{-1};
    amrex::Real m_type_forest;
    amrex::Real m_x_forest;
    amrex::Real m_y_forest;
    amrex::Real m_height_forest;
    amrex::Real m_diameter_forest;
    amrex::Real m_cd_forest;
    amrex::Real m_lai_forest;
    amrex::Real m_laimax_forest;

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real lm() const
    {
        amrex::Real treelaimax = 0.0;
        if (m_type_forest == 2) {
            const amrex::Real treeZm = m_laimax_forest * m_height_forest;
            const int n = 100;
            const amrex::Real int_0_fh = amr_wind::utils::trapz(
                0, m_height_forest - constants::TIGHT_TOL, n,
                [=](const amrex::Real x) noexcept {
                    const amrex::Real ratio =
                        (m_height_forest - treeZm) / (m_height_forest - x);
                    const auto exponent = (x < treeZm) ? 6.0 : 0.5;
                    return std::pow(ratio, exponent) *
                           std::exp(exponent * (1 - ratio));
                });
            treelaimax = m_lai_forest / int_0_fh;
        }
        return treelaimax;
    }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
    area_fraction(const amrex::Real z, const amrex::Real treelaimax) const
    {
        amrex::Real af = 0.0;
        if (m_type_forest == 1) {
            af = m_lai_forest / m_height_forest;
        } else if (m_type_forest == 2) {
            const auto treeZm = m_laimax_forest * m_height_forest;
            const auto ratio =
                (m_height_forest - treeZm) / (m_height_forest - z);
            const auto exponent = (z < treeZm) ? 6.0 : 0.5;
            af = treelaimax * std::pow(ratio, exponent) *
                 std::exp(exponent * (1 - ratio));
        }
        return af;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::RealBox real_bounding_box(
        const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& prob_lo,
        const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& prob_hi) const
    {
        const amrex::Real search_tol = 1.1;
        const amrex::Real search_radius = 0.5 * search_tol * m_diameter_forest;
        const auto x0 = m_x_forest - search_radius;
        const auto y0 = m_y_forest - search_radius;
        const auto z0 = prob_lo[2];
        const auto x1 = m_x_forest + search_radius;
        const auto y1 = m_y_forest + search_radius;
        const auto z1 = prob_hi[2];
        return {x0, y0, z0, x1, y1, z1};
    }

    amrex::Box bounding_box(const amrex::Geometry& geom) const
    {
        return utils::realbox_to_box(
            real_bounding_box(geom.ProbLoArray(), geom.ProbHiArray()), geom);
    }
};

/** Forestdrag Flow physics
 *  \ingroup physics
 */

class ForestDrag : public Physics::Register<ForestDrag>
{
public:
    static std::string identifier() { return "ForestDrag"; }

    explicit ForestDrag(CFDSim& sim);

    ~ForestDrag() override = default;

    void
    initialize_fields(int /*level*/, const amrex::Geometry& /*geom*/) override;

    void pre_init_actions() override {}

    void post_init_actions() override {}

    void post_regrid_actions() override;

    void pre_advance_work() override {}

    void post_advance_work() override {}

    amrex::Vector<Forest> read_forest(const int level) const;

private:
    CFDSim& m_sim;
    Field& m_forest_drag;
    Field& m_forest_id;
    std::string m_forest_file{"forest.amrwind"};
    bool m_verbose{false};

    /**
     * \brief Compute the global floating-point bounding box (amrex::RealBox)
     * that encloses all forests, using their physical extents.
     *
     * \param forests Vector of Forest objects to include in the bounding box.
     * \return        amrex::RealBox representing the global bounding box of all
     * forests.
     */
    static amrex::RealBox
    global_forest_realbox(const amrex::Vector<Forest>& forests)
    {
        amrex::Real xmin = std::numeric_limits<amrex::Real>::max();
        amrex::Real xmax = std::numeric_limits<amrex::Real>::lowest();
        amrex::Real ymin = std::numeric_limits<amrex::Real>::max();
        amrex::Real ymax = std::numeric_limits<amrex::Real>::lowest();
        amrex::Real hmax = std::numeric_limits<amrex::Real>::lowest();

        for (const auto& forest : forests) {
            const amrex::Real search_tol = 1.1;
            const amrex::Real search_radius =
                0.5 * search_tol * forest.m_diameter_forest;
            const amrex::Real xlo = forest.m_x_forest - search_radius;
            const amrex::Real xhi = forest.m_x_forest + search_radius;
            const amrex::Real ylo = forest.m_y_forest - search_radius;
            const amrex::Real yhi = forest.m_y_forest + search_radius;
            xmin = std::min(xmin, xlo);
            xmax = std::max(xmax, xhi);
            ymin = std::min(ymin, ylo);
            ymax = std::max(ymax, yhi);
            hmax = std::max(hmax, forest.m_height_forest);
        }

        return amrex::RealBox({xmin, ymin, 0.0}, {xmax, ymax, hmax});
    }

    /**
     * \brief Compute the maximum terrain height in a given box.
     *
     * \param arr  Array4 for terrain height.
     * \param box  AMReX box to search.
     * \return     Maximum terrain height in the box.
     */
    static amrex::Real max_terrain_in_box(
        const amrex::Array4<const amrex::Real>& arr, const amrex::Box& box)
    {
        amrex::Real max_terrain = std::numeric_limits<amrex::Real>::lowest();
        for (int k = box.smallEnd(2); k <= box.bigEnd(2); ++k) {
            for (int j = box.smallEnd(1); j <= box.bigEnd(1); ++j) {
                for (int i = box.smallEnd(0); i <= box.bigEnd(0); ++i) {
                    max_terrain = std::max(max_terrain, arr(i, j, k));
                }
            }
        }
        return max_terrain;
    }
};
} // namespace amr_wind::forestdrag

#endif
