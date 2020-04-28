#ifndef DENSITYSOURCE_H
#define DENSITYSOURCE_H

#include "Factory.H"
#include "FieldDescTypes.H"
#include "FieldUtils.H"
#include "FieldRepo.H"
#include "AMReX_MultiFab.H"

namespace amr_wind {

class CFDSim;

namespace pde {

class DensitySource : public Factory<DensitySource, const CFDSim&>
{
public:
    static const std::string base_identifier() { return "DensitySource"; }

    virtual ~DensitySource() = default;

    virtual void operator()(
        const int lev,
        const amrex::MFIter& mfi,
        const amrex::Box& bx,
        const FieldState fstate,
        const amrex::Array4<amrex::Real>& src_term) const = 0;
};

} // namespace pde
} // namespace amr_wind

#endif /* DENSITYSOURCE_H */
