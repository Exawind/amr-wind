#ifndef RainDrop_H
#define RainDrop_H

#include "amr-wind/core/Physics.H"
#include "amr-wind/core/Field.H"

/** Multiphase falling raindrop physics
 *  \ingroup mulitphase_physics
 *
 */

namespace amr_wind {

class RainDrop : public Physics::Register<RainDrop>
{
    static_assert(AMREX_SPACEDIM == 3, "RainDrop requires 3 dimensional mesh");

public:
    static std::string identifier() { return "RainDrop"; }

    explicit RainDrop(CFDSim& sim);

    ~RainDrop() override = default;

    //! Initialize the levelset and velocity fields for Rain Drop simulations

    void initialize_fields(int level, const amrex::Geometry& geom) override;

    void post_init_actions() override {}

    void post_regrid_actions() override {}

    void pre_advance_work() override {}

    void post_advance_work() override {}

private:
    Field& m_velocity;
    Field& m_levelset;

    //! Not const because calculation of VOF takes place during initialize
    CFDSim& m_sim;

    //! Initial Rain Drop location
    amrex::Vector<amrex::Real> m_loc{0.006, 0.006, 0.006};

    //! Initial Rain Drop velocity
    amrex::Vector<amrex::Real> m_vel{0.0, 0.0, 0.0};

    //! Raindrop radius value
    amrex::Real m_radius{0.003};
};

} // namespace amr_wind

#endif /* RainDrop_H */
