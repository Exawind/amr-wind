#include "AmrTestMesh.H"
#include "gtest/gtest.h"
#include "AMReX_MultiFab.H"

namespace amr_wind_tests {

TestLevelData::TestLevelData(
    const amrex::BoxArray& ba, const amrex::DistributionMapping& dm)
    : m_ba(ba), m_dm(dm), m_factory(new amrex::FArrayBoxFactory())
{}

amrex::MultiFab& TestLevelData::declare_field(
    const std::string& name, const int ncomp, const int num_ghost)
{
    auto found = m_fields.find(name);
    if (found != m_fields.end()) {
        amrex::Abort("Attempt to register existing field: " + name);
    }

    m_fields[name] = amrex::MultiFab();
    m_fields[name].define(
        m_ba, m_dm, ncomp, num_ghost, amrex::MFInfo(), *m_factory);

    return m_fields[name];
}

amrex::MultiFab& TestLevelData::get_field(const std::string& name)
{
    auto found = m_fields.find(name);
    if (found == m_fields.end()) {
        amrex::Abort("Cannot find field: " + name);
    }
    return m_fields[name];
}

amrex::Vector<amrex::MultiFab*> AmrTestMesh::declare_field(
    const std::string& name, const int ncomp, const int num_ghost)
{
    amrex::Vector<amrex::MultiFab*> vec(finest_level+1, nullptr);
    for (int lev=0; lev < (finest_level + 1); ++lev) {
        vec[lev] = &m_leveldata[lev]->declare_field(name, ncomp, num_ghost);
    }
    return vec;
}

amrex::Vector<amrex::MultiFab*> AmrTestMesh::get_field(const std::string& name)
{
    amrex::Vector<amrex::MultiFab*> vec(finest_level+1, nullptr);
    for (int lev=0; lev < (finest_level + 1); ++lev) {
        vec[lev] = &m_leveldata[lev]->get_field(name);
    }
    return vec;
}

AmrTestMesh::AmrTestMesh()
    : m_leveldata(max_level + 1),
      m_repo(*this)
{}

void AmrTestMesh::initialize_mesh(amrex::Real current_time)
{
    InitFromScratch(current_time);
}

void AmrTestMesh::MakeNewLevelFromScratch(
    int lev,
    amrex::Real time,
    const amrex::BoxArray& ba,
    const amrex::DistributionMapping& dm)
{
    SetBoxArray(lev, ba);
    SetDistributionMap(lev, dm);

    m_leveldata[lev].reset(
        new TestLevelData(boxArray(lev), DistributionMap(lev)));

    m_repo.make_new_level_from_scratch(lev, time, ba, dm);
}

void AmrTestMesh::MakeNewLevelFromCoarse(
    int lev,
    amrex::Real time,
    const amrex::BoxArray& ba,
    const amrex::DistributionMapping& dm)
{
    SetBoxArray(lev, ba);
    SetDistributionMap(lev, dm);

    // Concrete test fixtures must override and initialize field data if
    // necessary

    m_repo.make_new_level_from_coarse(lev, time, ba, dm);
}

void AmrTestMesh::RemakeLevel(
    int lev,
    amrex::Real time,
    const amrex::BoxArray& ba,
    const amrex::DistributionMapping& dm)
{
    SetBoxArray(lev, ba);
    SetDistributionMap(lev, dm);

    m_repo.remake_level(lev, time, ba, dm);
}

void AmrTestMesh::ClearLevel(int lev)
{
    m_repo.clear_level(lev);
}

void AmrTestMesh::ErrorEst(
    int /* lev */, amrex::TagBoxArray& /* tags */, amrex::Real /* time */, int /* ngrow */)
{
    amrex::Abort("Not implemented");
}

} // namespace amr_wind_tests
