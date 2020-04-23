#ifndef PHYSICS_H
#define PHYSICS_H

#include "Factory.H"
#include "CollMgr.H"
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

class Physics : public Factory<Physics, const CFDSim&>
{
public:
    static const std::string base_identifier() { return "Physics"; }

    virtual ~Physics() = default;

    //! Initialize the fields for a given level during init/regrid
    virtual void initialize_fields(int level, const amrex::Geometry& geom) = 0;

    //! Perform tasks necessary before advancing timestep
    virtual void pre_advance_work() = 0;
};

class PhysicsMgr : public CollMgr<PhysicsMgr, Physics>
{};

}

#endif /* PHYSICS_H */
