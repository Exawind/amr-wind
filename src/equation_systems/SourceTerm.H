#ifndef SOURCETERM_H
#define SOURCETERM_H

#include "Factory.H"
#include "FieldDescTypes.H"
#include "FieldUtils.H"
#include "FieldRepo.H"
#include "AMReX_MultiFab.H"

namespace amr_wind {
namespace pde {

class SourceTerm
{
public:
    virtual ~SourceTerm() = default;

    virtual void operator()(
        const int lev,
        const amrex::MFIter& mfi,
        const amrex::Box& bx,
        const FieldState fstate,
        const amrex::Array4<amrex::Real>& src_term
    ) const = 0;
};

}
}


#endif /* SOURCETERM_H */
