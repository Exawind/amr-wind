#ifndef SloshingTank_H
#define SloshingTank_H

#include "amr-wind/core/Physics.H"
#include "amr-wind/core/Field.H"

/** Multiphase Sloshing Tank physics
 *  \ingroup mulitphase_physics
 *
 */

namespace amr_wind {

class SloshingTank : public Physics::Register<SloshingTank>
{
public:
    static const std::string identifier() { return "SloshingTank"; }

    explicit SloshingTank(CFDSim& sim);

    virtual ~SloshingTank() = default;

    //! Initialize the levelset and velocity fields for Sloshing Tank
    //! simulations
    void initialize_fields(int level, const amrex::Geometry& geom) override;

    void post_init_actions() override {}

    void post_regrid_actions() override {}

    void pre_advance_work() override {}

    void post_advance_work() override {}

private:
    Field& m_velocity;
    Field& m_levelset;

    //! Initial free surface amplitude magnitude
    amrex::Real m_amplitude{0.05};

    //! Initial free surface wavelength
    amrex::Real m_wavelength{4.0};

    //! Initial zero-level free-surface water depth
    amrex::Real m_waterlevel{0.0};
};

} // namespace amr_wind
#endif /* SloshingTank_H */
