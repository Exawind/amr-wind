#ifndef PHYSICS_H
#define PHYSICS_H

#include "amr-wind/core/Factory.H"
#include "amr-wind/core/CollMgr.H"
#include "AMReX_MultiFab.H"
#include "AMReX_Geometry.H"

namespace amr_wind {

class CFDSim;

class PhysicsOld
{
public:
    virtual ~PhysicsOld() = default;

    virtual std::string physics_name() const = 0;

    //! Initialize the fields for the physics/problem
    virtual void initialize_fields(
        int level,
        const amrex::Geometry& geom) const = 0;

    //! Perform tasks necessary before advancing timestep
    virtual void pre_advance_work() = 0;

    //! Add momentum source terms
    virtual void add_momentum_sources(
        const amrex::Geometry& geom,
        const amrex::MultiFab& density,
        const amrex::MultiFab& velocity,
        const amrex::MultiFab& scalars,
        amrex::MultiFab& vel_forces) const = 0;
};

/**
 *  \defgroup physics Physics models
 *
 *  AMR-Wind representation of different physics
 *
 *  In AMR-Wind, different physics are implemented using the Physics class.
 *  Every simulation is a collection of several physics instances that are used
 *  to set initial conditions, perform certain actions at the start and
 *  completion of a timestep, or after a regrid.
 *
 *  Examples of pre-processing actions include reading time-varying inflow data
 *  from an input file based on user inputs.
 */

/** Abstract representation of different physics
 *  \ingroup physics
 *
 *  This class defines an abstract API that represents the notion of some
 *  physics that will be used during the simulation. The most common use-case
 *  for this class is to create initial conditions. However, it can also be used
 *  to perform actions at different stages within a timestep during a CFD
 *  simulation.
 */
class Physics : public Factory<Physics, CFDSim&>
{
public:
    static const std::string base_identifier() { return "Physics"; }

    virtual ~Physics() = default;

    //! Initialize the fields for a given level during init/regrid
    virtual void initialize_fields(int level, const amrex::Geometry& geom) = 0;

    //! Actions to perform post initialization
    virtual void post_init_actions() = 0;

    //! Actions to perform post regrid
    virtual void post_regrid_actions() = 0;

    //! Perform tasks necessary before advancing timestep
    virtual void pre_advance_work() = 0;

    //! Perform tasks necessary after advancing timestep
    virtual void post_advance_work() = 0;
};

/** A collection of \ref physics instances that are active during a simulation
 *
 *  \ingroup physics
 */
class PhysicsMgr : public CollMgr<PhysicsMgr, Physics>
{
public:
    //! Return a concrete physics instance
    template<typename T>
    T& get()
    {
        AMREX_ASSERT(contains(T::identifier()));
        return dynamic_cast<T&>(operator()(T::identifier()));
    }

    //! Return a concrete physics instance
    template<typename T>
    const T& get() const
    {
        AMREX_ASSERT(contains(T::identifier()));
        return dynamic_cast<const T&>(operator()(T::identifier()));
    }
};

}

#endif /* PHYSICS_H */
