#ifndef DamBreak_H
#define DamBreak_H

#include "amr-wind/core/Physics.H"
#include "amr-wind/core/Field.H"

/** Multiphase Dam break physics
 *  \ingroup mulitphase_physics
 *
 */

namespace amr_wind {

class DamBreak : public Physics::Register<DamBreak>
{
public:
    static const std::string identifier() { return "DamBreak"; }

    explicit DamBreak(CFDSim& sim);

    virtual ~DamBreak() = default;

    //! Initialize the levelset and velocity fields for Dam Break simulations

    void initialize_fields(int level, const amrex::Geometry& geom) override;

    void post_init_actions() override {}

    void post_regrid_actions() override {}

    void pre_advance_work() override {}

    void post_advance_work() override {}

private:
    Field& m_velocity;
    Field& m_levelset;

    //! Initial DamBreak location
    amrex::Vector<amrex::Real> m_loc{{0.0, 0.0, 0.0}};

    //! dam width value
    amrex::Real m_width{0.1};

    //! dam height value
    amrex::Real m_height{0.25};
};

} // namespace amr_wind

#endif /* DamBreak_H */
