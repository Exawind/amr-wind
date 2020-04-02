#ifndef FIELDREPO_H
#define FIELDREPO_H

#include <string>
#include <unordered_map>

#include "FieldDescTypes.H"
#include "FieldUtils.H"
#include "Field.H"
#include "ScratchField.H"

#include "AMReX_AmrCore.H"
#include "AMReX_MultiFab.H"

namespace amr_wind {

struct LevelDataHolder
{
    amrex::Vector<amrex::MultiFab> m_mfabs;
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
    ~FieldRepo() = default;

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

    inline Field& declare_cc_field(
        const std::string& name,
        const int ncomp = 1,
        const int ngrow = 0,
        const int nstates = 1)
    {
        return declare_field(name, ncomp, ngrow, nstates, FieldLoc::CELL);
    }

    inline Field& declare_nd_field(
        const std::string& name,
        const int ncomp = 1,
        const int ngrow = 0,
        const int nstates = 1)
    {
        return declare_field(name, ncomp, ngrow, nstates, FieldLoc::NODE);
    }

    inline Field& declare_xf_field(
        const std::string& name,
        const int ncomp = 1,
        const int ngrow = 0,
        const int nstates = 1)
    {
        return declare_field(name, ncomp, ngrow, nstates, FieldLoc::XFACE);
    }

    inline Field& declare_yf_field(
        const std::string& name,
        const int ncomp = 1,
        const int ngrow = 0,
        const int nstates = 1)
    {
        return declare_field(name, ncomp, ngrow, nstates, FieldLoc::YFACE);
    }

    inline Field& declare_zf_field(
        const std::string& name,
        const int ncomp = 1,
        const int ngrow = 0,
        const int nstates = 1)
    {
        return declare_field(name, ncomp, ngrow, nstates, FieldLoc::ZFACE);
    }

    inline amrex::Vector<Field*> declare_face_normal_field(
        const amrex::Vector<std::string>& names,
        const int ncomp = 1,
        const int ngrow = 0,
        const int nstates = 1)
    {
        AMREX_ASSERT(names.size() == AMREX_SPACEDIM);
        return {
            &declare_field(names[0], ncomp, ngrow, nstates, FieldLoc::XFACE),
            &declare_field(names[1], ncomp, ngrow, nstates, FieldLoc::YFACE),
            &declare_field(names[2], ncomp, ngrow, nstates, FieldLoc::ZFACE),
        };
    }

    Field& get_field(
        const std::string& name, const FieldState fstate = FieldState::New) const;

    std::unique_ptr<ScratchField> create_scratch_field(
        const std::string& name,
        const int ncomp = 1,
        const int ngrow = 0,
        const FieldLoc floc = FieldLoc::CELL) const;

    std::unique_ptr<ScratchField> create_scratch_field(
        const int ncomp = 1,
        const int ngrow = 0,
        const FieldLoc floc = FieldLoc::CELL) const;

    int num_active_levels() const noexcept { return m_mesh.finestLevel() + 1; }

    int num_fields() const noexcept { return m_field_vec.size(); }

    const amrex::Vector<std::unique_ptr<Field>>& fields() const noexcept { return m_field_vec; }

    inline const amrex::FabFactory<amrex::FArrayBox>&
    factory(int lev) const noexcept { return *m_factory[lev]; }

protected:
    inline amrex::MultiFab& get_multifab(const unsigned fid, const int lev) noexcept
    {
        BL_ASSERT(lev <= m_mesh.finestLevel());
        return m_leveldata[lev]->m_mfabs[fid];
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

    mutable amrex::Vector<std::unique_ptr<Field>> m_field_vec;
    std::unordered_map<std::string, size_t> m_fid_map;

    bool m_is_initialized{false};
};

} // namespace amr_wind

#endif /* FIELDREPO_H */
