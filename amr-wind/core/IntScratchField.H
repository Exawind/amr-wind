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

/** A temporary computational field
 *  \ingroup fields
 *
 *  An IntScratchField is similar to a Field in usage, but is temporary in nature.
 *  It is used as a scratch buffer to compute intermediate quantities. However,
 *  unlike fields these don't have multiple states, and cannot survive across a
 *  regrid. By default, FieldRepo returns a unique pointer to this instance and
 *  it is not safe to hold this pointer across timesteps.
 *
 *  At present, IntScratchField cannot be used for I/O and/or post-processing
 * utilities.
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

    /** Return a sub-view of the ScratchField instance
     */
    /*
    ViewField<IntScratchField> subview(const int scomp = 0, const int ncomp = 1)
    {
        return {*this, scomp, ncomp};
    }
    */

    void fillpatch(amrex::Real time) noexcept;
    void fillpatch(amrex::Real time, const amrex::IntVect& ng) noexcept;

    void fillpatch(
        int lev,
        amrex::Real time,
        amrex::iMultiFab& mfab,
        const amrex::IntVect& nghost) noexcept;

    void fillpatch(
        int lev,
        amrex::Real time,
        amrex::iMultiFab& mfab,
        const amrex::IntVect& nghost,
        amrex::Vector<amrex::BCRec>& bcrec) noexcept;

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

#endif /* SCRATCHFIELD_H */
