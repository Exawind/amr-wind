#ifndef INTFIELD_H
#define INTFIELD_H

#include <string>

#include "FieldDescTypes.H"

#include "AMReX_iMultiFab.H"

namespace amr_wind {

class FieldRepo;

class IntField
{
public:
    friend class FieldRepo;

    IntField(const IntField&) = delete;
    IntField& operator=(const IntField&) = delete;

    //! Name of the field
    inline const std::string& name() const { return m_name; }

    //! Unique integer ID for this field
    inline unsigned id() const { return m_id; }

    //! Number of components for this field
    inline int num_comp() const { return m_ncomp; }

    //! Number of ghost cells
    inline const amrex::IntVect& num_grow() const { return m_ngrow; }

    //! Location of the field
    inline FieldLoc field_location() const { return m_floc; }

    //! Reference to the FieldRepo that holds the fabs
    const FieldRepo& repo() const { return m_repo; }

    //! Access the FAB at a given level
    amrex::iMultiFab& operator()(int lev) noexcept;
    const amrex::iMultiFab& operator()(int lev) const noexcept;

    amrex::Vector<amrex::iMultiFab*> vec_ptrs() noexcept;

    amrex::Vector<const amrex::iMultiFab*> vec_const_ptrs() const noexcept;

    void setVal(int value) noexcept;

    void setVal(int value, int start_comp, int num_comp=1, int nghost=0) noexcept;

    void setVal(const amrex::Vector<int>& values, int nghost=0) noexcept;

protected:
    IntField(
        FieldRepo& repo,
        const std::string& name,
        const unsigned fid,
        const int ncomp = 1,
        const int ngrow = 1,
        const FieldLoc floc = FieldLoc::CELL);

    FieldRepo& m_repo;

    std::string m_name;

    const unsigned m_id;

    int m_ncomp;

    amrex::IntVect m_ngrow;

    FieldLoc m_floc;
};

} // namespace amr_wind

#endif /* INTFIELD_H */
