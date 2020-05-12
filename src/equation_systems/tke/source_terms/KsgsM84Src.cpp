#include "KsgsM84Src.H"
#include "CFDSim.H"
#include "TurbulenceModel.H"

namespace amr_wind {
namespace pde {
namespace tke {

KsgsM84Src::KsgsM84Src(const CFDSim& sim)
{
    AMREX_ALWAYS_ASSERT(sim.turbulence_model().model_name() == "OneEqKsgsM84");
    auto coeffs = sim.turbulence_model().model_coeffs();
    m_Ceps = coeffs["Ceps"];
}

KsgsM84Src::~KsgsM84Src() = default;

void KsgsM84Src::operator()(
    const int lev,
    const amrex::MFIter& mfi,
    const amrex::Box& bx,
    const FieldState fstate,
    const amrex::Array4<amrex::Real>& src_term) const
{}

}
}
}
