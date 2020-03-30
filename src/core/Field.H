#ifndef FIELD_H
#define FIELD_H

#include <string>
#include <memory>
#include <unordered_map>

#include "AMReX_MultiFab.H"

#include "FieldDescTypes.H"

namespace amr_wind {

class Field;
class FieldRepo;

struct FieldInfo
{
    FieldInfo(
        const std::string& basename,
        const int ncomp,
        const int ngrow,
        const int nstates,
        const FieldLoc floc)
        : m_basename(basename)
        , m_ncomp(ncomp)
        , m_ngrow(ngrow)
        , m_nstates(nstates)
        , m_floc(floc)
    {
    }

    ~FieldInfo() = default;

    std::string m_basename;
    int m_ncomp;
    amrex::IntVect m_ngrow;
    int m_nstates;
    FieldLoc m_floc;

    std::unordered_map<FieldState, Field*> m_states;
};

class Field
{
public:
    friend class FieldRepo;

    Field(const Field&) = delete;
    Field& operator=(const Field&) = delete;

    inline const std::string& name() const { return m_name; }
    inline int num_comp() const { return m_info->m_ncomp; }
    inline const amrex::IntVect& num_grow() const { return m_info->m_ngrow; }
    inline int num_states() const { return m_info->m_nstates; }
    inline FieldLoc field_location() const { return m_info->m_floc; }
    inline FieldState field_state() const { return m_state; }

    Field& state(const FieldState fstate);
    const Field& state(const FieldState fstate) const;

    amrex::MultiFab& operator()(int lev) noexcept;
    const amrex::MultiFab& operator()(int lev) const noexcept;

    amrex::Vector<amrex::MultiFab*> vec_ptrs() noexcept;

    amrex::Vector<const amrex::MultiFab*> vec_const_ptrs() const noexcept;

protected:
    Field(
        FieldRepo& repo,
        const std::string& name,
        const std::shared_ptr<FieldInfo>& finfo,
        FieldState state);

    //! Reference to the FieldRepository instance
    FieldRepo& m_repo;

    //! Name of the field
    std::string m_name;

    //! Common data for all field states
    std::shared_ptr<FieldInfo> m_info;

    //! State for this field
    FieldState m_state;
};

} // namespace amr_wind

#endif /* FIELD_H */
