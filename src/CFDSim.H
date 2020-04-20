#ifndef CFDSIM_H
#define CFDSIM_H

#include "AMReX_AmrCore.H"
#include "SimTime.H"
#include "FieldRepo.H"
#include "PDEBase.H"

namespace amr_wind {

/** Data structures for a CFD simulation
 *
 *  CFDSim is a thin wrapper that holds all the necessary objects for a CFD
 *  simulation. The key data members within this object are:
 *
 *  - mesh (amrex::AmrCore) The AMR mesh hierarchy data structure
 *  - time (SimTime)        The time object
 *  - repo (FieldRepo)      The field repository
 */
class CFDSim
{
public:
    CFDSim(amrex::AmrCore& mesh) : m_mesh(mesh), m_time(), m_repo(m_mesh) {}

    ~CFDSim() = default;

    //! Return the AMR mesh hierarchy
    amrex::AmrCore& mesh() { return m_mesh; }
    const amrex::AmrCore& mesh() const { return m_mesh; }

    //! Return simulation time control
    SimTime& time() { return m_time; }
    const SimTime& time() const { return m_time; }

    //! Return the field repository
    FieldRepo& repo() { return m_repo; }
    const FieldRepo& repo() const { return m_repo; }

private:
    amrex::AmrCore& m_mesh;

    SimTime m_time;

    FieldRepo m_repo;
};

}

#endif /* CFDSIM_H */
