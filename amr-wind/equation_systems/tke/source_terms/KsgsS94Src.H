#ifndef KSGSM84SRC_H
#define KSGSM84SRC_H

#include "tke/TKESource.H"

namespace amr_wind {
namespace pde {
namespace tke {

class KsgsS94Src : public TKESource::Register<KsgsS94Src>
{
public:
    static const std::string identifier() { return "KsgsS94Src"; }

    KsgsS94Src(const CFDSim&);

    virtual ~KsgsS94Src();

    virtual void operator()(
        const int lev,
        const amrex::MFIter& mfi,
        const amrex::Box& bx,
        const FieldState fstate,
        const amrex::Array4<amrex::Real>& src_term) const override;
};
}
}
}

#endif /* KSGSM84SRC_H */
