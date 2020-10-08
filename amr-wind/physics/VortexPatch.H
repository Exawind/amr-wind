#ifndef VortexPatch_H
#define VortexPatch_H

#include "amr-wind/core/Physics.H"
#include "amr-wind/physics/VortexPatchFieldInit.H"
#include "amr-wind/core/Field.H"

/** Vortex patch physics
 *  \ingroup physics
 *
 */

namespace amr_wind {

class VortexPatch : public Physics::Register<VortexPatch>
{
public:
    static const std::string identifier() { return "VortexPatch"; }

   VortexPatch(CFDSim& sim);

   virtual ~VortexPatch() = default;

   //! Initialize the levelset and velocity fields for Vortex Patch simulations

    void initialize_fields(
        int level,
        const amrex::Geometry& geom) override;

    void post_init_actions() override {}

    void post_regrid_actions() override {}

    void pre_advance_work() override {} 

    void post_advance_work() override;

private:
    const CFDSim& m_sim;

    Field& m_velocity;
    Field& m_density;

    //! Levelset field
    Field* m_levelset{nullptr};

    //! VortexPatch field initializer instance
    std::unique_ptr<VortexPatchFieldInit> m_field_init;

    //! prescribed velocity period
    amrex::Real m_TT{3.0};

};

}

#endif /* VortexPatch_H */
