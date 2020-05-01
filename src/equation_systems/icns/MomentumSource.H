#ifndef MOMENTUMSOURCE_H
#define MOMENTUMSOURCE_H

#include "Factory.H"
#include "FieldDescTypes.H"
#include "FieldUtils.H"
#include "FieldRepo.H"
#include "AMReX_MultiFab.H"

namespace amr_wind {

class CFDSim;

namespace pde {

class MomentumSource : public Factory<MomentumSource, const CFDSim&>
{
public:
    static const std::string base_identifier() { return "MomentumSource"; }

    virtual ~MomentumSource() = default;

    virtual void operator()(
        const int lev,
        const amrex::MFIter& mfi,
        const amrex::Box& bx,
        const FieldState fstate,
        const amrex::Array4<amrex::Real>& src_term) const = 0;
};

} // namespace pde
} // namespace amr_wind

#endif /* MOMENTUMSOURCE_H */
