#ifndef VIEWFIELD_H
#define VIEWFIELD_H

#include <string>

#include "amr-wind/core/FieldDescTypes.H"
#include "AMReX_MultiFab.H"
#include "AMReX_Vector.H"

namespace amr_wind {

class FieldRepo;

/** Subview of a Field or a ScratchField
 *  \ingroup fields
 *
 *  A ViewField instance holds MultiFab corresponding to a subview of a Field or
 *  a ScratchField. As the data is managed by the source Field or ScratchField,
 *  subview instances cannot be created directly but instead must be created
 *  using the Field::subview() and ScratchField::subview() methods.
 */
template <typename T>
class ViewField
{
public:
    friend T;

    ViewField(const ViewField&) = delete;
    ViewField& operator=(const ViewField&) = delete;

    ViewField(ViewField&&) = default;
    ViewField& operator=(ViewField&&) = default;

    inline const std::string name() const { return m_src.name() + "view"; }

    //! Number of components in this subview
    inline int num_comp() const { return m_ncomp; }

    inline const amrex::IntVect& num_grow() const { return m_src.num_grow(); }

    inline FieldLoc field_location() const { return m_src.field_location(); }

    //! Return the field data for a given level
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

    //! Return a reference to the field repository that created this field
    const FieldRepo& repo() const { return m_src.repo(); }

    //! Return a reference to the source field from which this view is derived
    T& source_field() { return m_src; }
    const T& source_field() const { return m_src; }

private:
    ViewField(T& src, const int scomp = 0, const int ncomp = 1);

    //! The base field
    T& m_src;

    //! Index of the starting component
    int m_scomp;

    //! Number of sub-components in this view
    int m_ncomp;

    //! MultiFabs holding the alias views
    amrex::Vector<amrex::MultiFab> m_data;
};

} // namespace amr_wind

#endif /* VIEWFIELD_H */
