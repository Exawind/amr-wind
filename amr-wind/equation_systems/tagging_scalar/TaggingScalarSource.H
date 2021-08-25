#ifndef TAGGINGSCALARSOURCE_H
#define TAGGINGSCALARSOURCE_H

#include "amr-wind/core/Factory.H"
#include "amr-wind/core/FieldDescTypes.H"
#include "amr-wind/core/FieldUtils.H"
#include "amr-wind/core/FieldRepo.H"
#include "AMReX_MultiFab.H"

namespace amr_wind {

class CFDSim;

namespace pde {

/** Representation of a tagging scalar source term
 *  \ingroup temp_src
 *
 *  All tagging scalar source terms must inherit from this class
 */
class TaggingScalarSource : public Factory<TaggingScalarSource, const CFDSim&>
{
public:
    static const std::string base_identifier() { return "TaggingScalarSource"; }

    virtual ~TaggingScalarSource() = default;

    virtual void operator()(
        const int lev,
        const amrex::MFIter& mfi,
        const amrex::Box& bx,
        const FieldState fstate,
        const amrex::Array4<amrex::Real>& src_term) const = 0;
};

} // namespace pde
} // namespace amr_wind

#endif /* TAGGINGSCALARSOURCE_H */
