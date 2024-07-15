#ifndef VariableTemperature_H
#define VariableTemperature_H

#include "amr-wind/core/Physics.H"
#include "amr-wind/core/Field.H"
#include "amr-wind/CFDSim.H"

namespace amr_wind::VariableTemperature {

namespace {} // namespace

/** VariableTemperature Flow physics
 *  \ingroup physics
 Provides a way to read non-uniform surface temperature into AMR - Wind
 */

class VariableTemperature : public Physics::Register<VariableTemperature>
{
public:
    static std::string identifier() { return "VariableTemperature"; }

    explicit VariableTemperature(CFDSim& sim);

    ~VariableTemperature() override = default;

    void
    initialize_fields(int /*level*/, const amrex::Geometry& /*geom*/) override
    {}

    void pre_init_actions() override;

    void post_init_actions() override;

    void post_regrid_actions() override {}

    void pre_advance_work() override;

    void post_advance_work() override {}

private:
    const amr_wind::SimTime& m_time;
    CFDSim& m_sim;
    const FieldRepo& m_repo;
    const amrex::AmrCore& m_mesh;
    Field& m_velocity;
    // Surface Temperature Field
    Field& m_surf_temp;
    amrex::Real m_current_time{100};
    amrex::Real m_separation{100};
    amrex::Vector<amrex::Real> m_xloc;
    amrex::Vector<amrex::Real> m_yloc;
    amrex::Vector<amrex::Real> m_Tloc;
    amrex::Gpu::DeviceVector<amrex::Real> device_xloc;
    amrex::Gpu::DeviceVector<amrex::Real> device_yloc;
    amrex::Gpu::DeviceVector<amrex::Real> device_Tloc;
};
} // namespace amr_wind::VariableTemperature

#endif
