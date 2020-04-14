#ifndef ADVOP_MOL_H
#define ADVOP_MOL_H

#include <type_traits>

#include "MOL.H"
#include "SchemeTraits.H"
#include "PDETraits.H"
#include "PDEOps.H"

#include "AMReX_Gpu.H"
#include "AMReX_ParmParse.H"

namespace amr_wind {
namespace pde {

template <typename PDE>
struct AdvectionOp<
    PDE,
    fvm::MOL,
    typename std::enable_if<std::is_base_of<ScalarTransport, PDE>::value>::type>
{
    AdvectionOp(PDEFields& fields_in) : fields(fields_in)
    {}

    void operator()(const FieldState , const amrex::Real )
    {
        static_assert(PDE::ndim == 1, "Invalid number of components for scalar");
    }

    PDEFields& fields;
};

} // namespace pde
} // namespace amr_wind



#endif /* ADVOP_MOL_H */
