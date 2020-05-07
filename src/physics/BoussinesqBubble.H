#ifndef BoussinesqBubble_H
#define BoussinesqBubble_H

#include "Physics.H"
#include "BoussinesqBubbleFieldInit.H"
#include "Field.H"

namespace amr_wind {

/** Boussinesq Bubble physics
 */
class BoussinesqBubble : public Physics::Register<BoussinesqBubble>
{
public:
    static const std::string identifier() { return "BoussinesqBubble"; }

    BoussinesqBubble(const CFDSim& sim);

    virtual ~BoussinesqBubble() = default;

    //! Initialize the temperature and velocity fields for BoussinesqBubble simulations
    void initialize_fields(
        int level,
        const amrex::Geometry& geom) override;

    void post_init_actions() override;

    void post_regrid_actions() override {}

    void pre_advance_work() override {}

    void post_advance_work() override {}

private:
    Field& m_velocity;
    Field& m_density;
    Field& m_temperature;

    //! BoussinesqBubble field initializer instance
    std::unique_ptr<BoussinesqBubbleFieldInit> m_field_init;
};

}

#endif /* BoussinesqBubble_H */
