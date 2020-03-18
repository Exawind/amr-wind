#ifndef PHYSICS_H
#define PHYSICS_H

#include "AMReX_MultiFab.H"
#include "AMReX_Geometry.H"

namespace amr_wind {

class Physics
{
public:
    virtual ~Physics() = default;

    virtual std::string physics_name() const = 0;

    //! Initialize the fields for the physics/problem
    virtual void initialize_fields(
        const amrex::Geometry& geom,
        amrex::MultiFab& density,
        amrex::MultiFab& velocity,
        amrex::MultiFab& pressure,
        amrex::MultiFab& scalars) const = 0;
};

}

#endif /* PHYSICS_H */
