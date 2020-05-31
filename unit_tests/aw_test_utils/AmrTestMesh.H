#ifndef AMRTESTMESH_H
#define AMRTESTMESH_H

#include <memory>
#include <unordered_map>

#include "AMReX_AmrCore.H"
#include "amr-wind/CFDSim.H"
#include "amr-wind/core/FieldRepo.H"
#include "amr-wind/core/Field.H"

namespace amr_wind_tests {

/** Base class for tests that require a valid Mesh
 *
 *  This class specializes amrex::AmrCore and provides a barebones interface
 *  that can be used within unit test fixtures.
 */
class AmrTestMesh : public amrex::AmrCore
{
public:
    AmrTestMesh();

    virtual ~AmrTestMesh() = default;

    //! Create the initial AMR hierarchy
    virtual void initialize_mesh(amrex::Real current_time);

    //! Return the total number of existing levels in AMR hierarchy
    int num_levels() const { return finest_level + 1; }

    amr_wind::CFDSim& sim() { return m_sim; }

    amr_wind::FieldRepo& field_repo() { return m_repo; }

protected:
    virtual void MakeNewLevelFromScratch(
        int lev, amrex::Real time, const amrex::BoxArray& ba,
        const amrex::DistributionMapping& dm) override;

    virtual void MakeNewLevelFromCoarse(
        int lev, amrex::Real time, const amrex::BoxArray& ba,
        const amrex::DistributionMapping& dm) override;

    virtual void RemakeLevel(
        int lev, amrex::Real time, const amrex::BoxArray& ba,
        const amrex::DistributionMapping& dm) override;

    virtual void ClearLevel(int lev) override;

    virtual void
    ErrorEst(int lev, amrex::TagBoxArray& tags, amrex::Real time, int ngrow)
        override;

    amr_wind::CFDSim m_sim;
    amr_wind::FieldRepo& m_repo;
};

}

#endif /* AMRTESTMESH_H */
