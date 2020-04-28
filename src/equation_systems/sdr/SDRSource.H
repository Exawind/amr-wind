#ifndef SDRSOURCE_H
#define SDRSOURCE_H

#include "Factory.H"
#include "FieldDescTypes.H"
#include "FieldUtils.H"
#include "FieldRepo.H"
#include "AMReX_MultiFab.H"

namespace amr_wind {

class CFDSim;

namespace pde {

class SDRSource : public Factory<SDRSource, const CFDSim&>
{
public:
    static const std::string base_identifier() { return "SDRSource"; }

    virtual ~SDRSource() = default;

    virtual void operator()(
        const int lev,
        const amrex::MFIter& mfi,
        const amrex::Box& bx,
        const FieldState fstate,
        const amrex::Array4<amrex::Real>& src_term) const = 0;
};

} // namespace pde
} // namespace amr_wind

#endif /* SDRSOURCE_H */
