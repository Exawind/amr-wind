#ifndef ABL_H
#define ABL_H

#include <memory>

#include "Physics.H"
#include "ABLFieldInit.H"
#include "ABLWallFunction.H"
#include "ABLForcing.H"
#include "BoussinesqBuoyancy.H"
#include "CoriolisForcing.H"
#include "SimTime.H"
#include "PlaneAveraging.H"

class incflo;

namespace amr_wind {

namespace pde {
namespace icns {
class ABLForcing;
}
}

/** Atmospheric Boundary Layer physics
 */
class ABLOld : public PhysicsOld
{
public:
    ABLOld(const SimTime& time, incflo* incflo_in);

    ~ABLOld() = default;

    //! Unique name indentifying this physics set for debugging
    std::string physics_name() const override { return "ABL"; }

    //! Initialize the temperature and velocity fields for ABL simulations
    void initialize_fields(
        int level,
        const amrex::Geometry& geom) const override;

    void pre_advance_work() override;

    //! Add momentum source terms
    void add_momentum_sources(
        const amrex::Geometry& geom,
        const amrex::MultiFab& density,
        const amrex::MultiFab& velocity,
        const amrex::MultiFab& scalars,
        amrex::MultiFab& vel_forces) const override;

    const ABLFieldInit& abl_initializer() const { return *m_field_init; }

private:
    const SimTime& m_time;

    // TODO: deprecate this
    incflo* m_incflo;

    PlaneAveraging m_pa;

    //! ABL field initializer instance
    std::unique_ptr<ABLFieldInit> m_field_init;

    //! Wall function instance
    std::unique_ptr<ABLWallFunction> m_abl_wall_func;

    //! ABL forcing function instance
    std::unique_ptr<ABLForcingOld> m_abl_forcing;

    //! Boussinesq buoyancy forcing term
    std::unique_ptr<BoussinesqBuoyancyOld> m_boussinesq;

    //! Coriolis forcing instance
    std::unique_ptr<CoriolisForcingOld> m_coriolis;

    //! True if Boussinesq buoyancy term is active
    bool m_has_boussinesq{true};

    //! True if Coriolis forcing term is active
    bool m_has_coriolis{true};

    //! True if driving pressure gradient term is activate
    bool m_has_driving_dpdx{true};
};

class ABL : public Physics::Register<ABL>
{
public:
    static const std::string identifier() { return "ABL"; }

    ABL(const CFDSim& sim);

    virtual ~ABL();

    void initialize_fields(
        int level,
        const amrex::Geometry& geom) override;

    void pre_advance_work() override;

    void register_forcing_term(pde::icns::ABLForcing* forcing) const
    {
        m_abl_forcing = forcing;
    }

private:
    const CFDSim& m_sim;

    Field& m_velocity;
    Field& m_density;
    Field& m_temperature;

    PlaneAveraging m_pa;

    //! ABL field initializer instance
    std::unique_ptr<ABLFieldInit> m_field_init;

    //! Reference to ABL forcing term if present
    mutable pde::icns::ABLForcing* m_abl_forcing{nullptr};
};

}

#endif /* ABL_H */
