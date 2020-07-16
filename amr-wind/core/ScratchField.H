#ifndef SCRATCHFIELD_H
#define SCRATCHFIELD_H

#include <string>

#include "amr-wind/core/FieldDescTypes.H"
#include "AMReX_MultiFab.H"

namespace amr_wind {

class FieldRepo;

/** A temporary computational field
 *  \ingroup fields
 *
 *  A ScratchField is similar to a Field in usage, but is temporary in nature.
 *  It is used as a scratch buffer to compute intermediate quantities. However,
 *  unlike fields these don't have multiple states, and cannot survive across a
 *  regrid. By default, FieldRepo returns a unique pointer to this instance and
 *  it is not safe to hold this pointer across timesteps.
 *
 *  At present, ScratchField cannot be used for I/O and/or post-processing utilities.
 */
class ScratchField
{
public:
    friend class FieldRepo;

    ScratchField(const ScratchField&) = delete;
    ScratchField& operator=(const ScratchField&) = delete;

    inline const std::string& name() const { return m_name; }

    //! Number of components for this field
    inline int num_comp() const { return m_ncomp; }

    //! Ghost cells
    inline const amrex::IntVect& num_grow() const { return m_ngrow; }

    inline FieldLoc field_location() const { return m_floc; }

    amrex::MultiFab& operator()(int lev) { return m_data[lev]; }
    const amrex::MultiFab& operator()(int lev) const { return m_data[lev]; }

    amrex::Vector<amrex::MultiFab*> vec_ptrs() noexcept
    {
        return amrex::GetVecOfPtrs(m_data);
    }

    amrex::Vector<const amrex::MultiFab*> vec_const_ptrs() const noexcept
    {
        return amrex::GetVecOfConstPtrs(m_data);
    }

    const FieldRepo& repo() const { return m_repo; }

protected:
    ScratchField(
        const FieldRepo& repo,
        const std::string& name,
        const int ncomp = 1,
        const int ngrow = 1,
        const FieldLoc floc = FieldLoc::CELL)
        : m_repo(repo), m_name(name), m_ncomp(ncomp), m_ngrow(ngrow), m_floc(floc)
    {}

    const FieldRepo& m_repo;
    std::string m_name;
    int m_ncomp;
    amrex::IntVect m_ngrow;
    FieldLoc m_floc;

    amrex::Vector<amrex::MultiFab> m_data;
};

}

#endif /* SCRATCHFIELD_H */
