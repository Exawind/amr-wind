#ifndef MetMastMarking_H
#define MetMastMarking_H

#include "amr-wind/core/Physics.H"
#include "amr-wind/core/Field.H"
#include "amr-wind/CFDSim.H"

namespace amr_wind::MetMastMarking {

/** MetMastMarking Flow physics
 *  \ingroup physics
 */

class MetMastMarking : public Physics::Register<MetMastMarking>
{
public:
    static std::string identifier() { return "MetMastMarking"; }

    explicit MetMastMarking(CFDSim& sim);

    ~MetMastMarking() override = default;

    void initialize_fields(int level, const amrex::Geometry& geom) override;

    void pre_init_actions() override {}

    void post_init_actions() override {}

    void post_regrid_actions() override;

    void pre_advance_work() override {}

    void post_advance_work() override {}

private:
    CFDSim& m_sim;
    const FieldRepo& m_repo;
    const amrex::AmrCore& m_mesh;

    //! Blanking field for met mast
    Field& m_metmast_blank;

    amrex::Vector<amrex::Real> m_metmast_x;
    amrex::Vector<amrex::Real> m_metmast_y;
    amrex::Vector<amrex::Real> m_metmast_z;
    amrex::Vector<amrex::Real> m_metmast_horizontal_radius;
    amrex::Vector<amrex::Real> m_metmast_vertical_radius;
    amrex::Gpu::DeviceVector<amrex::Real> m_metmast_x_d;
    amrex::Gpu::DeviceVector<amrex::Real> m_metmast_y_d;
    amrex::Gpu::DeviceVector<amrex::Real> m_metmast_z_d;
    amrex::Gpu::DeviceVector<amrex::Real> m_metmast_horizontal_radius_d;
    amrex::Gpu::DeviceVector<amrex::Real> m_metmast_vertical_radius_d;
    //! Met mast file
    std::string m_metmast_file{"metmast.amrwind"};
};
} // namespace amr_wind::MetMastMarking

#endif
