#ifndef REFINEMENTCRITERIA_H
#define REFINEMENTCRITERIA_H

#include "AMReX_AmrCore.H"
#include "AMReX_TagBox.H"

namespace amr_wind {

//! Abstract interface for refinement criteria
class RefinementCriteria
{
public:
    RefinementCriteria() = default;

    virtual ~RefinementCriteria() = default;

    //! Perform initialization tasks (e.g., read from file etc.). Called before
    //! mesh is initialized
    virtual void initialize(const amrex::AmrCore&) = 0;

    /** Perform tagging operation based on refinement criteria
     *
     *  Interface is the same as amrex::AmrCore::ErrorEst
     */
    virtual void operator()(
        int level,
        amrex::TagBoxArray& tags,
        amrex::Real time,
        int ngrow) = 0;
};

}

#endif /* REFINEMENTCRITERIA_H */
