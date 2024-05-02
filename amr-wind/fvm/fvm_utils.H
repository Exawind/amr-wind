#ifndef FVM_UTILS_H
#define FVM_UTILS_H

#include "amr-wind/core/FieldRepo.H"
#include "amr-wind/fvm/stencils.H"

#include "AMReX_MFIter.H"

namespace amr_wind::fvm::impl {

/** Apply a finite volume operator for a given field
 */
template <typename FvmOp, typename FType>
inline void apply(const FvmOp& fvmop, const FType& fld)
{
    namespace stencil = amr_wind::fvm::stencil;
    const int nlevels = fld.repo().num_active_levels();
    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& domain = fld.repo().mesh().Geom(lev).Domain();
        const auto& mfab = fld(lev);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        for (amrex::MFIter mfi(mfab, amrex::TilingIfNotGPU()); mfi.isValid();
             ++mfi) {

            fvmop.template apply<stencil::StencilInterior>(lev, mfi);

            // Check if the box touches any of the physical domain boundaries.
            // If not, short circuit the rest of the code and move on to the
            // next multifab.
            if (domain.strictly_contains(mfi.tilebox())) {
                continue;
            }

            // faces
            fvmop.template apply<stencil::StencilILO>(lev, mfi);
            fvmop.template apply<stencil::StencilJLO>(lev, mfi);
            fvmop.template apply<stencil::StencilKLO>(lev, mfi);
            fvmop.template apply<stencil::StencilIHI>(lev, mfi);
            fvmop.template apply<stencil::StencilJHI>(lev, mfi);
            fvmop.template apply<stencil::StencilKHI>(lev, mfi);

            // edges
            fvmop.template apply<stencil::StencilIHI_JLO>(lev, mfi);
            fvmop.template apply<stencil::StencilIHI_JHI>(lev, mfi);

            fvmop.template apply<stencil::StencilIHI_KLO>(lev, mfi);
            fvmop.template apply<stencil::StencilIHI_KHI>(lev, mfi);

            fvmop.template apply<stencil::StencilJHI_KLO>(lev, mfi);
            fvmop.template apply<stencil::StencilJHI_KHI>(lev, mfi);

            fvmop.template apply<stencil::StencilILO_JLO>(lev, mfi);
            fvmop.template apply<stencil::StencilILO_JHI>(lev, mfi);

            fvmop.template apply<stencil::StencilILO_KLO>(lev, mfi);
            fvmop.template apply<stencil::StencilILO_KHI>(lev, mfi);

            fvmop.template apply<stencil::StencilJLO_KLO>(lev, mfi);
            fvmop.template apply<stencil::StencilJLO_KHI>(lev, mfi);

            // corners
            fvmop.template apply<stencil::StencilILO_JLO_KLO>(lev, mfi);
            fvmop.template apply<stencil::StencilILO_JLO_KHI>(lev, mfi);
            fvmop.template apply<stencil::StencilILO_JHI_KLO>(lev, mfi);
            fvmop.template apply<stencil::StencilILO_JHI_KHI>(lev, mfi);
            fvmop.template apply<stencil::StencilIHI_JLO_KLO>(lev, mfi);
            fvmop.template apply<stencil::StencilIHI_JLO_KHI>(lev, mfi);
            fvmop.template apply<stencil::StencilIHI_JHI_KLO>(lev, mfi);
            fvmop.template apply<stencil::StencilIHI_JHI_KHI>(lev, mfi);
        }
    }
}

} // namespace amr_wind::fvm::impl

#endif /* FVM_UTILS_H */
