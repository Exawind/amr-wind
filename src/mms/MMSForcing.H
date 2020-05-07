#ifndef MMSFORCING_H
#define MMSFORCING_H

#include "FieldRepo.H"
#include "icns/MomentumSource.H"

namespace amr_wind {
namespace pde {
namespace icns {
namespace mms {

class MMSForcing : public MomentumSource::Register<MMSForcing>
{
public:
    static const std::string identifier() { return "MMSForcing"; }

    MMSForcing(const CFDSim&);

    virtual ~MMSForcing();

    virtual void operator()(
        const int lev,
        const amrex::MFIter& mfi,
        const amrex::Box& bx,
        const FieldState fstate,
        const amrex::Array4<amrex::Real>& src_term) const override;

private:
    const Field& m_mms_vel_source;
};
} // namespace mms
} // namespace icns
} // namespace pde
} // namespace amr_wind

#endif /* MMSFORCING_H */
