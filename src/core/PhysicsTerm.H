#ifndef PHYSICSTERM_H
#define PHYSICSTERM_H

#include "AMReX_MultiFab.H"
#include "AMReX_Geometry.H"
#include "LevelData.H"

namespace amr_wind {

class PhysicsTerm
{
public:
    virtual ~PhysicsTerm() = default;

    virtual std::string physics_name() const = 0;

    //! Initialize fields before starting the simulation
    virtual void initialize_fields(const int lev) = 0;

    //! Perform tasks necessary before advancing timestep
    virtual void pre_advance_work() = 0;

    //! Add momentum source terms
    virtual void add_momentum_sources(
        const int lev,
        amrex::MultiFab& vel_forces) const = 0;
};

}

#endif /* PHYSICSTERM_H */
