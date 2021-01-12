#ifndef VortexPatch_H
#define VortexPatch_H

#include "amr-wind/core/Physics.H"
#include "amr-wind/core/Field.H"

/** Multiphase vortex patch physics
 *  \ingroup mulitphase_physics
 *
 */

namespace amr_wind {

class VortexPatch : public Physics::Register<VortexPatch>
{
    static_assert(
        AMREX_SPACEDIM == 3, "VortexPatch requires 3 dimensional mesh");

public:
    static const std::string identifier() { return "VortexPatch"; }

    explicit VortexPatch(CFDSim& sim);

    virtual ~VortexPatch() = default;

    //! Initialize the levelset and velocity fields for Vortex Patch simulations

    void initialize_fields(int level, const amrex::Geometry& geom) override;

    void post_init_actions() override {}

    void post_regrid_actions() override {}

    void pre_advance_work() override;

    void post_advance_work() override;

private:
    const CFDSim& m_sim;

    Field& m_velocity;
    Field& m_levelset;
    Field& m_density;

    //! Initial VortexPatch location
    amrex::Vector<amrex::Real> m_loc{{0.35, 0.35, 0.35}};

    //! vortex patch radius value
    amrex::Real m_radius{0.15};

    //! prescribed velocity period
    amrex::Real m_TT{3.0};
};

} // namespace amr_wind

#endif /* VortexPatch_H */
