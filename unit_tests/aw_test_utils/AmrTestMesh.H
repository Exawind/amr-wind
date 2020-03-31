#ifndef AMRTESTMESH_H
#define AMRTESTMESH_H

#include <memory>
#include <unordered_map>

#include "AMReX_AmrCore.H"
#include "FieldRepo.H"
#include "Field.H"

namespace amr_wind_tests {

/** Manager for handling data at a given AMR level
 */
class TestLevelData
{
public:
    using FabFactoryType = amrex::FArrayBoxFactory;

    TestLevelData(
        const amrex::BoxArray& ba,
        const amrex::DistributionMapping& dm);

    virtual ~TestLevelData() = default;

    /** Create field data at this level
     *
     *  @param name Unique name to refer to this field
     *  @param ncomp Number of components (default 1)
     *  @param num_ghost Number of ghost cells (default 0)
     *
     *  @return MultiFab for the newly declared field
     */
    virtual amrex::MultiFab&
    declare_field(const std::string& name, const int ncomp = 1, const int num_ghost = 0);

    /** Return reference to a previously created field
     *
     *  @param name Unique string identifier for the field
     */
    virtual amrex::MultiFab& get_field(const std::string& name);

protected:
    const amrex::BoxArray& m_ba;
    const amrex::DistributionMapping& m_dm;

    std::unique_ptr<amrex::FabFactory<amrex::FArrayBox>> m_factory;

    std::unordered_map<std::string, amrex::MultiFab> m_fields;
};

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

    /** Create field data on on all levels of the AmrMesh
     *
     *  @param name Unique name to refer to this field
     *  @param ncomp Number of components (default 1)
     *  @param num_ghost Number of ghost cells (default 0)
     *
     *  @return MultiFab for the newly declared field
     */
    amrex::Vector<amrex::MultiFab*> declare_field(
        const std::string& name, const int ncomp = 1, const int num_ghost = 0);

    /** Return a vector of existing field on all levels
     *
     *  @param name Unique string identifier for the field
     */
    amrex::Vector<amrex::MultiFab*> get_field(const std::string& name);

    //! Return the total number of existing levels in AMR hierarchy
    int num_levels() const { return finest_level + 1; }

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

    amrex::Vector<std::unique_ptr<TestLevelData>> m_leveldata;

    amr_wind::FieldRepo m_repo;
};

}

#endif /* AMRTESTMESH_H */
