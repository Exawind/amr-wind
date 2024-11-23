#ifndef ForestDrag_H
#define ForestDrag_H

#include "amr-wind/core/Physics.H"
#include "amr-wind/core/Field.H"
#include "amr-wind/CFDSim.H"

namespace amr_wind::Forestdrag {

namespace {} // namespace

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
    struct Forest
    {
        amrex::Real type_forest;
        amrex::Real x_forest;
        amrex::Real y_forest;
        amrex::Real height_forest;
        amrex::Real diameter_forest;
        amrex::Real cd_forest;
        amrex::Real lai_forest;
        amrex::Real laimax_forest;
    };
    amrex::Vector<Forest> m_forests;
    amrex::Gpu::DeviceVector<Forest> m_d_forests;
};
} // namespace amr_wind::Forestdrag

#endif
