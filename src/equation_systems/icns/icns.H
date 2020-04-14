#ifndef ICNS_H
#define ICNS_H

#include "PDETraits.H"
#include "SchemeTraits.H"
#include "PDEHelpers.H"
#include "PDE.H"
#include "icns/MomentumSource.H"

namespace amr_wind {
namespace pde {

struct ICNS : VectorTransport
{
    using MLDiffOp = amrex::MLTensorOp;
    using SrcTerm = MomentumSource;

    static std::string pde_name() { return "ICNS"; }
    static std::string var_name() { return "velocity"; }
};

} // namespace pde
} // namespace amr_wind

#endif /* ICNS_H */
