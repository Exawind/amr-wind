#ifndef SOURCETERM_H
#define SOURCETERM_H

#include "amr-wind/core/Factory.H"
#include "amr-wind/core/FieldDescTypes.H"
#include "amr-wind/core/FieldUtils.H"
#include "amr-wind/core/FieldRepo.H"
#include "AMReX_MultiFab.H"

namespace amr_wind {
class CFDSim;

namespace pde {

class SourceTerm : public Factory<SourceTerm, const CFDSim&>
{
public:
    static const std::string base_identifier() { return "SourceTerm"; }

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
