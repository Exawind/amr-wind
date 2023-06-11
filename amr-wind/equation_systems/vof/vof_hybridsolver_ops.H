#ifndef VOF_HYBRIDSOLVER_OPS_H_
#define VOF_HYBRIDSOLVER_OPS_H_

#include <AMReX_FArrayBox.H>
#include "amr-wind/core/FieldRepo.H"
#include <cmath>

namespace amr_wind::multiphase {

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real sharpen_kernel(
    int i,
    int j,
    int k,
    amrex::Array4<amrex::Real const> const& volfrac) noexcept
{
    return std::max(
        0.0,
        std::min(
            1.0,
            0.5 + (volfrac(i, j, k) < 0.5 ? -1.0 : 1.0) *
                      std::pow(std::abs(volfrac(i, j, k) - 0.5), 1.0 / 3.0)));
}

static void sharpen_acquired_vof(
    const int nlevels, amr_wind::IntField& f_iblank, amr_wind::Field& f_vof)
{
    // Sharpen data from nalu-wind (in iblank regions)
    for (int lev = 0; lev < nlevels; ++lev) {
        auto& iblank = f_iblank(lev);
        auto& vof = f_vof(lev);

        for (amrex::MFIter mfi(iblank); mfi.isValid(); ++mfi) {
            const auto& gbx = mfi.growntilebox();
            const amrex::Array4<const int>& native_flag =
                iblank.const_array(mfi);
            const amrex::Array4<amrex::Real>& volfrac = vof.array(mfi);
            amrex::ParallelFor(
                gbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    // In iblanked regions, sharpen VOF and limit it
                    volfrac(i, j, k) = (native_flag(i, j, k) > 0)
                                           ? volfrac(i, j, k)
                                           : sharpen_kernel(i, j, k, volfrac);
                });
        }
    }
}

} // namespace amr_wind::multiphase

#endif // VOF_HYBRIDSOLVER_OPS.H
