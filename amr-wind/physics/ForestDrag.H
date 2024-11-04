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
    amrex::Vector<amrex::Real> m_type_forest;
    amrex::Vector<amrex::Real> m_x_forest;
    amrex::Vector<amrex::Real> m_y_forest;
    amrex::Vector<amrex::Real> m_ht_forest;
    amrex::Vector<amrex::Real> m_dia_forest;
    amrex::Vector<amrex::Real> m_cd_forest;
    amrex::Vector<amrex::Real> m_lai_forest;
    amrex::Vector<amrex::Real> m_lm_forest;
};
} // namespace amr_wind::Forestdrag

#endif
