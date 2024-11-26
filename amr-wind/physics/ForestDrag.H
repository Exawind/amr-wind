#ifndef ForestDrag_H
#define ForestDrag_H

#include "amr-wind/core/Physics.H"
#include "amr-wind/core/Field.H"
#include "amr-wind/CFDSim.H"
#include "amr-wind/utilities/index_operations.H"

namespace amr_wind::forestdrag {

namespace {
struct Forest
{
    amrex::Real m_type_forest;
    amrex::Real m_x_forest;
    amrex::Real m_y_forest;
    amrex::Real m_height_forest;
    amrex::Real m_diameter_forest;
    amrex::Real m_cd_forest;
    amrex::Real m_lai_forest;
    amrex::Real m_laimax_forest;

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real calc_lm() const
    {
        amrex::Real treelaimax = 0.0;
        if (m_type_forest == 2) {
            amrex::Real ztree = 0;
            amrex::Real expFun = 0;
            amrex::Real ratio = 0;
            const amrex::Real treeZm = m_laimax_forest * m_height_forest;
            const amrex::Real dz = m_height_forest / 100;
            while (ztree <= m_height_forest) {
                ratio = (m_height_forest - treeZm) / (m_height_forest - ztree);
                if (ztree < treeZm) {
                    expFun = expFun +
                             std::pow(ratio, 6.0) * std::exp(6 * (1 - ratio));
                } else {
                    expFun = expFun +
                             std::pow(ratio, 0.5) * std::exp(0.5 * (1 - ratio));
                }
                ztree = ztree + dz;
            }
            treelaimax = m_lai_forest / (expFun * dz);
        }
        return treelaimax;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
    calc_af(const amrex::Real z, const amrex::Real treelaimax) const
    {
        if (m_type_forest == 1) {
            return m_lai_forest / m_height_forest;
        } else if (m_type_forest == 2) {
            const auto treeZm = m_laimax_forest * m_height_forest;
            const auto ratio =
                (m_height_forest - treeZm) / (m_height_forest - z);
            const auto exponent = (z < treeZm) ? 6.0 : 0.5;
            return treelaimax * std::pow(ratio, exponent) *
                   std::exp(exponent * (1 - ratio));
        }
        return 0.0;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::RealBox real_bounding_box(
        const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& prob_lo) const
    {
        const auto x1 = m_x_forest - 0.5 * m_diameter_forest;
        const auto y1 = m_y_forest - 0.5 * m_diameter_forest;
        const auto z1 = prob_lo[2];
        const auto x2 = m_x_forest + 0.5 * m_diameter_forest;
        const auto y2 = m_y_forest + 0.5 * m_diameter_forest;
        const auto z2 = m_height_forest;
        return {x1, y1, z1, x2, y2, z2};
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Box
    bounding_box(const amrex::Geometry& geom) const
    {
        return utils::realbox_to_box(
            real_bounding_box(geom.ProbLoArray()), geom);
    }
};
} // namespace

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

private:
    CFDSim& m_sim;
    Field& m_forest_drag;
    Field& m_forest_blank;
    amrex::Vector<Forest> m_forests;
    amrex::Gpu::DeviceVector<Forest> m_d_forests;
};
} // namespace amr_wind::forestdrag

#endif
