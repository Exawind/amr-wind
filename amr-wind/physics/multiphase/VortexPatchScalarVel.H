#ifndef VortexPatchScalarVel_H
#define VortexPatchScalarVel_H

#include "amr-wind/core/Physics.H"
#include "amr-wind/core/Field.H"

/** Multiphase vortex patch scalar velocity physics
 *  \ingroup multiphase_physics
 *
 */

//! This test case is intended as an evaluation of the momentum advection
//! scheme. The cell-centered velocity field is treated as a scalar, initialized
//! as a smooth function crossing the interface. The advection velocities at the
//! faces are prescribed as functions of time and space, and do not depend on
//! the changing cell-centered velocity field.

namespace amr_wind {

class VortexPatchScalarVel : public Physics::Register<VortexPatchScalarVel>
{
    static_assert(
        AMREX_SPACEDIM == 3,
        "VortexPatchScalarVel requires 3 dimensional mesh");

public:
    static std::string identifier() { return "VortexPatchScalarVel"; }

    explicit VortexPatchScalarVel(CFDSim& sim);

    ~VortexPatchScalarVel() override = default;

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

    //! Initial VortexPatchScalarVel location
    amrex::Vector<amrex::Real> m_loc{{0.35, 0.35, 0.35}};

    //! vortex patch radius value
    amrex::Real m_radius{0.15};

    //! smoothing factor
    amrex::Real m_sfactor{1.0};

    //! prescribed velocity period
    amrex::Real m_TT{3.0};
};

} // namespace amr_wind

#endif /* VortexPatchScalarVel_H */
