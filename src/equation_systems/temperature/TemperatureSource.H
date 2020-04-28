#ifndef TEMPERATURESOURCE_H
#define TEMPERATURESOURCE_H

#include "Factory.H"
#include "FieldDescTypes.H"
#include "FieldUtils.H"
#include "FieldRepo.H"
#include "AMReX_MultiFab.H"

namespace amr_wind {

class CFDSim;

namespace pde {

class TemperatureSource : public Factory<TemperatureSource, const CFDSim&>
{
public:
    static const std::string base_identifier() { return "TemperatureSource"; }

    virtual ~TemperatureSource() = default;

    virtual void operator()(
        const int lev,
        const amrex::MFIter& mfi,
        const amrex::Box& bx,
        const FieldState fstate,
        const amrex::Array4<amrex::Real>& src_term) const = 0;
};

} // namespace pde
} // namespace amr_wind


#endif /* TEMPERATURESOURCE_H */
