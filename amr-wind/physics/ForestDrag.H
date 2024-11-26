#ifndef ForestDrag_H
#define ForestDrag_H

#include "amr-wind/core/Physics.H"
#include "amr-wind/core/Field.H"
#include "amr-wind/CFDSim.H"

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

    amrex::Real calc_lm() const;
    amrex::Real
    calc_af(const amrex::Real z, const amrex::Real treelaimax) const;
};

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real Forest::calc_lm() const
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
                expFun =
                    expFun + std::pow(ratio, 6.0) * std::exp(6 * (1 - ratio));
            } else {
                expFun =
                    expFun + std::pow(ratio, 0.5) * std::exp(0.5 * (1 - ratio));
            }
            ztree = ztree + dz;
        }
        treelaimax = m_lai_forest / (expFun * dz);
    }
    return treelaimax;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
Forest::calc_af(const amrex::Real z, const amrex::Real treelaimax) const
{
    amrex::Real af = 0.0;
    if (m_type_forest == 1) {
        af = m_lai_forest / m_height_forest;
    } else if (m_type_forest == 2) {
        const amrex::Real treeZm = m_laimax_forest * m_height_forest;
        const amrex::Real ratio =
            (m_height_forest - treeZm) / (m_height_forest - z);
        if (z < treeZm) {
            af = treelaimax * std::pow(ratio, 6.0) * std::exp(6 * (1 - ratio));
        } else if (z <= m_height_forest) {
            af =
                treelaimax * std::pow(ratio, 0.5) * std::exp(0.5 * (1 - ratio));
        }
    }
    return af;
}
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
    initialize_fields(int /*level*/, const amrex::Geometry& /*geom*/) override
    {}

    void pre_init_actions() override;

    void post_init_actions() override;

    void post_regrid_actions() override {}

    void pre_advance_work() override {}

    void post_advance_work() override {}

    int return_blank_value(int i, int j, int k);

private:
    CFDSim& m_sim;
    const FieldRepo& m_repo;
    const amrex::AmrCore& m_mesh;
    Field& m_velocity;
    //! Forest Drag Force Term
    Field& m_forest_drag;
    Field& m_forest_blank;
    amrex::Vector<Forest> m_forests;
    amrex::Gpu::DeviceVector<Forest> m_d_forests;
};
} // namespace amr_wind::forestdrag

#endif
