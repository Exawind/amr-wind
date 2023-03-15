#ifndef INTSCRATCHFIELD_H
#define INTSCRATCHFIELD_H

#include <string>
#include <utility>

#include "amr-wind/core/FieldDescTypes.H"
#include "amr-wind/core/ViewField.H"
#include "AMReX_iMultiFab.H"
#include "AMReX_Vector.H"
#include "AMReX_PhysBCFunct.H"

namespace amr_wind {

class FieldRepo;

/** Integer version of ScratchField that works with iMultiFab
 */
class IntScratchField
{
public:
    friend class FieldRepo;

    IntScratchField(const IntScratchField&) = delete;
    IntScratchField& operator=(const IntScratchField&) = delete;

    //! Name if available for this scratch field
    inline const std::string& name() const { return m_name; }

    //! Number of components for this field
    inline int num_comp() const { return m_ncomp; }

    //! Ghost cells
    inline const amrex::IntVect& num_grow() const { return m_ngrow; }

    //! Cell, node, face where the field is stored
    inline FieldLoc field_location() const { return m_floc; }

    //! Return the field data for a given level
    amrex::iMultiFab& operator()(int lev) { return m_data[lev]; }
    const amrex::iMultiFab& operator()(int lev) const { return m_data[lev]; }

    amrex::Vector<amrex::iMultiFab*> vec_ptrs() noexcept
    {
        return amrex::GetVecOfPtrs(m_data);
    }

    amrex::Vector<const amrex::iMultiFab*> vec_const_ptrs() const noexcept
    {
        return amrex::GetVecOfConstPtrs(m_data);
    }

    //! Return a reference to the field repository that created this field
    const FieldRepo& repo() const { return m_repo; }

    void setVal(int value) noexcept;

protected:
    IntScratchField(
        const FieldRepo& repo,
        std::string name,
        const int ncomp = 1,
        const int ngrow = 1,
        const FieldLoc floc = FieldLoc::CELL)
        : m_repo(repo)
        , m_name(std::move(name))
        , m_ncomp(ncomp)
        , m_ngrow(ngrow)
        , m_floc(floc)
    {}

    const FieldRepo& m_repo;
    std::string m_name;
    int m_ncomp;
    amrex::IntVect m_ngrow;
    FieldLoc m_floc;

    amrex::Vector<amrex::iMultiFab> m_data;
};

} // namespace amr_wind

#endif /* INTSCRATCHFIELD_H */
