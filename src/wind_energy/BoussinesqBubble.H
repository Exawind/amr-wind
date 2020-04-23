#ifndef BoussinesqBubble_H
#define BoussinesqBubble_H

#include <memory>

#include "Physics.H"
#include "BoussinesqBubbleFieldInit.H"
#include "BoussinesqBuoyancy.H"
#include "SimTime.H"

class incflo;

namespace amr_wind {

/** Boussinesq Bubble physics
 */
class BoussinesqBubbleOld : public PhysicsOld
{
public:
    BoussinesqBubbleOld(incflo* incflo_in);

    ~BoussinesqBubbleOld() = default;

    //! Unique name indentifying this physics set for debugging
    std::string physics_name() const override { return "BoussinesqBubble"; }

    //! Initialize the temperature and velocity fields for BoussinesqBubble simulations
    void initialize_fields(
        int level,
        const amrex::Geometry& geom) const override;

    void pre_advance_work() override {};

    //! Add momentum source terms
    void add_momentum_sources(
        const amrex::Geometry& geom,
        const amrex::MultiFab& density,
        const amrex::MultiFab& velocity,
        const amrex::MultiFab& scalars,
        amrex::MultiFab& vel_forces) const override;

    const BoussinesqBubbleFieldInit& BoussinesqBubble_initializer() const { return *m_field_init; }

private:

    // TODO: deprecate this
    incflo* m_incflo;

    //! BoussinesqBubble field initializer instance
    std::unique_ptr<BoussinesqBubbleFieldInit> m_field_init;

    //! Boussinesq buoyancy forcing term
    std::unique_ptr<BoussinesqBuoyancyOld> m_boussinesq;

    //! True if Boussinesq buoyancy term is active
    bool m_has_boussinesq{true};

};

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

    void pre_advance_work() override {}

private:
    Field& m_velocity;
    Field& m_density;
    Field& m_temperature;

    //! BoussinesqBubble field initializer instance
    std::unique_ptr<BoussinesqBubbleFieldInit> m_field_init;
};

}

#endif /* BoussinesqBubble_H */
