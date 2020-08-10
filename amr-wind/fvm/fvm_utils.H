#ifndef FVM_UTILS_H
#define FVM_UTILS_H

#include "amr-wind/core/FieldRepo.H"
#include "amr-wind/fvm/stencils.H"

#include "AMReX_MFIter.H"

namespace amr_wind {
namespace fvm {
namespace impl {

template <typename FvmOp, typename FType>
inline void apply(const FvmOp& fvmop, const FType& fld)
{
    namespace stencil = amr_wind::fvm::stencil;
    const int nlevels = fld.repo().num_active_levels();
    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& mfab = fld(lev);
        for (amrex::MFIter mfi(mfab); mfi.isValid(); ++mfi) {
            fvmop.template apply<stencil::StencilInterior>(lev, mfi);
            fvmop.template apply<stencil::StencilILO>(lev, mfi);
            fvmop.template apply<stencil::StencilJLO>(lev, mfi);
            fvmop.template apply<stencil::StencilKLO>(lev, mfi);
            fvmop.template apply<stencil::StencilIHI>(lev, mfi);
            fvmop.template apply<stencil::StencilJHI>(lev, mfi);
            fvmop.template apply<stencil::StencilKHI>(lev, mfi);
        }
    }
}

}
}
}

#endif /* FVM_UTILS_H */
