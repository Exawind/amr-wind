#ifndef FIELDREPO_H
#define FIELDREPO_H

#include <string>
#include <unordered_map>

#include "FieldDescTypes.H"
#include "FieldUtils.H"
#include "Field.H"

#include "AMReX_AmrCore.H"
#include "AMReX_MultiFab.H"

namespace amr_wind {

struct LevelDataHolder
{
    std::unordered_map<std::string, std::unique_ptr<amrex::MultiFab>> m_data;
};

class FieldRepo
{
public:
    friend class Field;

    using FieldMap = std::unordered_map<std::string, std::unique_ptr<Field>>;

    FieldRepo(const amrex::AmrCore& mesh)
        : m_mesh(mesh)
        , m_leveldata(mesh.maxLevel() + 1)
        , m_factory(mesh.maxLevel() + 1)
    {}

    FieldRepo(const FieldRepo&) = delete;
    FieldRepo& operator=(const FieldRepo&) = delete;

    void make_new_level_from_scratch(
        int lev, amrex::Real time,
        const amrex::BoxArray& ba,
        const amrex::DistributionMapping& dm);

    void make_new_level_from_coarse(
        int lev, amrex::Real time,
        const amrex::BoxArray& ba,
        const amrex::DistributionMapping& dm);

    void remake_level(
        int lev, amrex::Real time,
        const amrex::BoxArray& ba,
        const amrex::DistributionMapping& dm);

    void clear_level(int lev);

    Field& declare_field(
        const std::string& name,
        const int ncomp = 1,
        const int ngrow = 0,
        const int nstates = 1,
        const FieldLoc floc = FieldLoc::CELL);

    Field& get_field(
        const std::string& name, const FieldState fstate = FieldState::New);

    int num_active_levels() const noexcept { return m_mesh.finestLevel() + 1; }

    const FieldMap& fields() const { return m_fields; }

protected:
    inline amrex::MultiFab* get_multifab(const std::string& name, int lev) noexcept
    {
        BL_ASSERT(lev <= m_mesh.finestLevel());
        return m_leveldata[lev]->m_data[name].get();
    }

    void allocate_field_data(
        int lev,
        Field& field,
        LevelDataHolder& fdata,
        const amrex::FabFactory<amrex::FArrayBox>& factory);

    void allocate_field_data(Field& field);

    void allocate_field_data(
        const amrex::BoxArray& ba,
        const amrex::DistributionMapping& dm,
        LevelDataHolder& fdata,
        const amrex::FabFactory<amrex::FArrayBox>& factory);

    const amrex::AmrCore& m_mesh;

    amrex::Vector<std::unique_ptr<LevelDataHolder>> m_leveldata;
    amrex::Vector<std::unique_ptr<amrex::FabFactory<amrex::FArrayBox>>> m_factory;

    FieldMap m_fields;

    bool m_is_initialized{false};
};

} // namespace amr_wind

#endif /* FIELDREPO_H */
