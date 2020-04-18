#ifndef PDETRAITS_H
#define PDETRAITS_H

#include <string>
#include "AMReX_MLTensorOp.H"

#include "SourceTerm.H"

namespace amr_wind {
namespace pde {

struct VectorTransport
{
    using MLDiffOp = amrex::MLABecLaplacian;
    static constexpr int ndim = AMREX_SPACEDIM;

    static constexpr bool multiply_rho = true;
    static constexpr bool has_diffusion = true;
};

struct ScalarTransport
{
    using MLDiffOp = amrex::MLABecLaplacian;
    using SrcTerm = SourceTerm;

    static constexpr int ndim = 1;
    static constexpr bool multiply_rho = true;
    static constexpr bool has_diffusion = true;
};

} // namespace pde
} // namespace amr_wind

#endif /* PDETRAITS_H */
