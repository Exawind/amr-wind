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

    static constexpr bool need_nph_state = true;
    static constexpr bool multiply_rho = true;
    static constexpr bool has_diffusion = true;
};

struct ScalarTransport
{
    //! AMReX linear operator used to solve the diffusion equation
    using MLDiffOp = amrex::MLABecLaplacian;
    //! Base class of the source term used to create specific instances
    using SrcTerm = SourceTerm;

    //! Number of components for the PDE
    static constexpr int ndim = 1;

    //! Do the PDE terms have to be multiplied by density
    static constexpr bool multiply_rho = true;

    //! Does this PDE have a diffusion term
    static constexpr bool has_diffusion = true;

    //! Default value to initialize the BC array
    static constexpr amrex::Real default_bc_value = 0.0;

    //! Does this scalar need an NPH state
    static constexpr bool need_nph_state = true;
};

} // namespace pde
} // namespace amr_wind

#endif /* PDETRAITS_H */
