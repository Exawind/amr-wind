#ifndef VOF_HYBRIDSOLVER_OPS_H_
#define VOF_HYBRIDSOLVER_OPS_H_

#include <AMReX_FArrayBox.H>
#include "amr-wind/core/FieldRepo.H"
#include <cmath>

namespace amr_wind::multiphase {

static void replace_masked_vof(
    const int nlevels,
    amr_wind::IntField& f_iblank,
    amr_wind::Field& f_vof,
    amr_wind::Field& f_vof_new)
{
    // Sharpen data from nalu-wind (in iblank regions)
    for (int lev = 0; lev < nlevels; ++lev) {
        auto& iblank = f_iblank(lev);
        auto& vof = f_vof(lev);
        const auto& vof_new = f_vof_new(lev);

        for (amrex::MFIter mfi(iblank); mfi.isValid(); ++mfi) {
            const auto& gbx = mfi.growntilebox();
            const amrex::Array4<const int>& native_flag =
                iblank.const_array(mfi);
            const amrex::Array4<amrex::Real>& volfrac = vof.array(mfi);
            const amrex::Array4<const amrex::Real>& vfmasked =
                vof_new.const_array(mfi);
            amrex::ParallelFor(
                gbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    // In iblanked regions, sharpen VOF and limit it
                    volfrac(i, j, k) = (native_flag(i, j, k) > 0)
                                           ? volfrac(i, j, k)
                                           : vfmasked(i, j, k);
                });
        }
    }
}

} // namespace amr_wind::multiphase

#endif // VOF_HYBRIDSOLVER_OPS.H
