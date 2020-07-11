#ifndef REFINEMENTCRITERIA_H
#define REFINEMENTCRITERIA_H

#include "AMReX_AmrCore.H"
#include "AMReX_TagBox.H"

/**
 *  \defgroup amr_utils Mesh refinement
 *  Mesh refinement
 *
 *  This group contains utilities used for static or adaptive refinement during
 *  a simulation.
 *
 *  \ingroup utilities
 */

namespace amr_wind {

/** Abstract interface for tagging cells for refinement
 *  \ingroup amr_utils
 *
 *  This class provides an API that can be used by concrete implementations to
 *  tag cells that must be refined based on a pre-defined criteria. The criteria
 *  can be a heuristic determined from the solution (e.g., vorticity magntiude
 *  or gradients), or some user defined criteria (e.g., static, nested refinements).
 */
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
