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
class BoussinesqBubble : public Physics
{
public:
    BoussinesqBubble(incflo* incflo_in);

    ~BoussinesqBubble() = default;

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
    std::unique_ptr<BoussinesqBuoyancy> m_boussinesq;

    //! True if Boussinesq buoyancy term is active
    bool m_has_boussinesq{true};

};

}

#endif /* BoussinesqBubble_H */
