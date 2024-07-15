#ifndef VariableRoughness_H
#define VariableRoughness_H

#include "amr-wind/core/Physics.H"
#include "amr-wind/core/Field.H"
#include "amr-wind/CFDSim.H"
#include "AMReX_ParmParse.H"
namespace amr_wind::VariableRoughness {

namespace {} // namespace

/** VariableRoughness Flow physics
 *  \ingroup physics
 Provides a way to include non uniform roughness length for ABL
 */

class VariableRoughness : public Physics::Register<VariableRoughness>
{
public:
    static std::string identifier() { return "VariableRoughness"; }

    explicit VariableRoughness(CFDSim& sim);

    ~VariableRoughness() override = default;

    void
    initialize_fields(int /*level*/, const amrex::Geometry& /*geom*/) override
    {}

    void pre_init_actions() override;

    void post_init_actions() override;

    void post_regrid_actions() override {}

    void pre_advance_work() override {}

    void post_advance_work() override {}

private:
    CFDSim& m_sim;
    const FieldRepo& m_repo;
    const amrex::AmrCore& m_mesh;
    Field& m_velocity;
    // Roughness Field
    Field& m_terrainz0;
    // Reading the Roughness Coordinates from file
    amrex::Vector<amrex::Real> m_xrough;
    amrex::Vector<amrex::Real> m_yrough;
    amrex::Vector<amrex::Real> m_z0rough;
};
} // namespace amr_wind::VariableRoughness

#endif
