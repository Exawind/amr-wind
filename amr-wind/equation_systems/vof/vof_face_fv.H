#ifndef VOF_FACE_FV_H_
#define VOF_FACE_FV_H_

#include <AMReX_MultiFabUtil.H>
#include "amr-wind/equation_systems/vof/volume_fractions.H"

namespace amr_wind {
namespace multiphase {

void set_density_face_via_vof(
    const amrex::Array<amrex::MultiFab*, AMREX_SPACEDIM>& fc,
    amrex::MultiFab& vof,
    amrex::MultiFab& VH,
    const amrex::Real rho1,
    const amrex::Real rho2)
{
    const amrex::Real tiny = 1e-12;

    for (amrex::MFIter mfi(vof); mfi.isValid(); ++mfi) {
        const auto& vbx = mfi.validbox();
        const auto& bx = mfi.growntilebox(1);
        const auto& F = vof.array(mfi);
        const auto& vof_half = VH.array(mfi);
        const auto& rho_xf = fc[0]->array(mfi);
        const auto& rho_yf = fc[1]->array(mfi);
        const auto& rho_zf = fc[2]->array(mfi);
        // Loop gets sub-cell phase volumes (cell loop)
        amrex::ParallelFor(
            bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                // Local VOF value and array for octant VOF
                const amrex::Real volfrac = F(i, j, k);
                amrex::GpuArray<amrex::Real, 8> subVOF;
                // Determine if cell is single-phase
                if (std::abs(volfrac - 1.0) <= tiny) {
                    // All liquid
                    for (int n = 0; n < 8; ++n) {
                        subVOF[n] = 1.0;
                    }
                } else if (std::abs(volfrac) <= tiny) {
                    // All gas
                    for (int n = 0; n < 8; ++n) {
                        subVOF[n] = 0.0;
                    }
                } else {
                    amrex::Real mx, my, mz, alpha;
                    // Get interface plane, according to current VOF field
                    fit_plane(i, j, k, F, mx, my, mz, alpha);
                    // Get subcell volumes through cutting
                    subVOF = cut_volume_oct(mx, my, mz, alpha);
                }
                // Store subcell volumes as half cell VOF
                // Left half (-x)
                vof_half(i, j, k, 0) =
                    0.25 * (subVOF[0] + subVOF[2] + subVOF[4] + subVOF[6]);
                // Right half (+x)
                vof_half(i, j, k, 1) =
                    0.25 * (subVOF[1] + subVOF[3] + subVOF[5] + subVOF[7]);
                // Bottom half (-y)
                vof_half(i, j, k, 2) =
                    0.25 * (subVOF[0] + subVOF[1] + subVOF[4] + subVOF[5]);
                // Top half (+y)
                vof_half(i, j, k, 3) =
                    0.25 * (subVOF[2] + subVOF[3] + subVOF[6] + subVOF[7]);
                // Back half (-z)
                vof_half(i, j, k, 4) =
                    0.25 * (subVOF[0] + subVOF[1] + subVOF[2] + subVOF[3]);
                // Front half (+z)
                vof_half(i, j, k, 5) =
                    0.25 * (subVOF[4] + subVOF[5] + subVOF[6] + subVOF[7]);
            });
        // Loop gets density of control volumes at faces (face loop)
        amrex::ParallelFor(
            bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                // At current cell, use half-cell VOF field to get left,
                // bottom, and top face densities. Use dedicated boxes to assign
                // to faces, avoid going out of bounds.
                amrex::Box const& xbx = amrex::surroundingNodes(vbx, 0);
                amrex::Box const& ybx = amrex::surroundingNodes(vbx, 1);
                amrex::Box const& zbx = amrex::surroundingNodes(vbx, 2);

                // NOTE: This step assumes constant density in each phase
                // and constant cell size on the current level.
                amrex::Real vof_fc;
                if (xbx.contains(i, j, k)) {
                    vof_fc =
                        0.5 * (vof_half(i - 1, j, k, 1) + vof_half(i, j, k, 0));
                    rho_xf(i, j, k) = vof_fc * rho1 + (1.0 - vof_fc) * rho2;
                }
                if (ybx.contains(i, j, k)) {
                    vof_fc =
                        0.5 * (vof_half(i, j - 1, k, 3) + vof_half(i, j, k, 2));
                    rho_yf(i, j, k) = vof_fc * rho1 + (1.0 - vof_fc) * rho2;
                }
                if (zbx.contains(i, j, k)) {
                    vof_fc =
                        0.5 * (vof_half(i, j, k - 1, 5) + vof_half(i, j, k, 4));
                    rho_zf(i, j, k) = vof_fc * rho1 + (1.0 - vof_fc) * rho2;
                }
            });
    }
}
} // namespace multiphase
} // namespace amr_wind

#endif
