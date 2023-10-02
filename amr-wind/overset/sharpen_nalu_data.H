#ifndef SHARPEN_NALU_DATA_H_
#define SHARPEN_NALU_DATA_H_

#include "amr-wind/CFDSim.H"
#include "amr-wind/physics/multiphase/MultiPhase.H"
#include "amr-wind/equation_systems/vof/volume_fractions.H"
#include "amr-wind/utilities/IOManager.H"

namespace amr_wind {
namespace overset {

// Approximate signed distance function
amrex::Real asdf(const amrex::Real a_vof, const amrex::Real i_th)
{
    // function of local vof value and interface thickness
    return (i_th * log((a_vof + 1e-12) / (1. - a_vof + 1e-12)));
}

void process_vof(amrex::MultiFab& mf_vof)
{
    constexpr amrex::Real tiny = 1e-6;
    for (amrex::MFIter mfi(mf_vof); mfi.isValid(); ++mfi) {
        const auto& gbx = mfi.growntilebox();
        const amrex::Array4<amrex::Real>& vof = mf_vof.array(mfi);
        // Populate approximate signed distance function
        amrex::ParallelFor(
            gbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                // This conditional is just for the sake of debugging
                if (vof(i, j, k) < 1.0 && vof(i, j, k) > 0.0) {
                    vof(i, j, k) =
                        vof(i, j, k) < tiny
                            ? 0.0
                            : (vof(i, j, k) > 1. - tiny ? 1. : vof(i, j, k));
                }
            });
    }
}

void populate_psi(
    amrex::MultiFab& mf_psi, amrex::MultiFab& mf_vof, const amrex::Real i_th)
{
    for (amrex::MFIter mfi(mf_psi); mfi.isValid(); ++mfi) {
        const auto& gbx = mfi.growntilebox();
        const amrex::Array4<amrex::Real>& psi = mf_psi.array(mfi);
        const amrex::Array4<const amrex::Real>& vof = mf_vof.const_array(mfi);
        // Populate approximate signed distance function
        amrex::ParallelFor(
            gbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                psi(i, j, k) = asdf(vof(i, j, k), i_th);
            });
    }
}

void populate_normal_vector(
    amrex::MultiFab& mf_normvec,
    amrex::MultiFab& mf_vof,
    amrex::iMultiFab& mf_iblank)
{
    for (amrex::MFIter mfi(mf_vof); mfi.isValid(); ++mfi) {
        const auto& gbxm1 = grow(mfi.growntilebox(), -1);
        const amrex::Array4<amrex::Real>& normvec = mf_normvec.array(mfi);
        const amrex::Array4<const amrex::Real>& vof = mf_vof.const_array(mfi);
        const amrex::Array4<const int>& iblank = mf_iblank.const_array(mfi);
        // Calculate gradients in each direction with centered diff
        // (Should this be replaced with youngs_fd_normal?)
        amrex::ParallelFor(
            gbxm1, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                // Neumann condition across nalu bdy
                int ibdy = 0;
                int jbdy = 0;
                int kbdy = 0;
                if (iblank(i, j, k) != iblank(i - 1, j, k)) {
                    ibdy = -1;
                }
                if (iblank(i, j, k) != iblank(i, j - 1, k)) {
                    jbdy = -1;
                }
                if (iblank(i, j, k) != iblank(i, j, k - 1)) {
                    kbdy = -1;
                }
                // no cell should be isolated such that -1 and 1 are needed
                if (iblank(i, j, k) != iblank(i + 1, j, k)) {
                    ibdy = +1;
                }
                if (iblank(i, j, k) != iblank(i, j + 1, k)) {
                    jbdy = +1;
                }
                if (iblank(i, j, k) != iblank(i, j, k + 1)) {
                    kbdy = +1;
                }
                // Calculate normal
                amrex::Real mx, my, mz, mmag;
                multiphase::youngs_fd_normal_neumann(
                    i, j, k, ibdy, jbdy, kbdy, vof, mx, my, mz);
                // Normalize normal
                mmag = std::sqrt(mx * mx + my * my + mz * mz + 1e-20);
                // Save normal
                normvec(i, j, k, 0) = mx / mmag;
                normvec(i, j, k, 1) = my / mmag;
                normvec(i, j, k, 2) = mz / mmag;
            });
    }
}

void populate_grad_psi(
    amrex::MultiFab& mf_gpsi,
    amrex::MultiFab& mf_vof,
    amrex::iMultiFab& mf_iblank,
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dx,
    const amrex::Real i_th)
{
    for (amrex::MFIter mfi(mf_gpsi); mfi.isValid(); ++mfi) {
        const auto& gbx = mfi.growntilebox();
        const auto& gbxm1 = grow(gbx, -1);
        const amrex::Array4<amrex::Real>& gpsi = mf_gpsi.array(mfi);
        const amrex::Array4<const amrex::Real>& vof = mf_vof.const_array(mfi);
        const amrex::Array4<const int>& iblank = mf_iblank.const_array(mfi);
        // Populate approximate signed distance function
        amrex::ParallelFor(
            gbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                gpsi(i, j, k, 3) = asdf(vof(i, j, k), i_th);
            });
        // Calculate gradients in each direction with centered diff
        // (Should this be replaced with youngs_fd_normal?)
        amrex::ParallelFor(
            gbxm1, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                // Neumann condition across nalu bdy
                if (iblank(i, j, k) == -1) {
                    const amrex::Real center = gpsi(i, j, k, 3);
                    const amrex::Real left = iblank(i - 1, j, k) == -1
                                                 ? gpsi(i - 1, j, k, 3)
                                                 : center;
                    const amrex::Real right = iblank(i + 1, j, k) == -1
                                                  ? gpsi(i + 1, j, k, 3)
                                                  : center;
                    const amrex::Real back = iblank(i, j - 1, k) == -1
                                                 ? gpsi(i, j - 1, k, 3)
                                                 : center;
                    const amrex::Real front = iblank(i, j + 1, k) == -1
                                                  ? gpsi(i, j + 1, k, 3)
                                                  : center;
                    const amrex::Real bottom = iblank(i, j, k - 1) == -1
                                                   ? gpsi(i, j, k - 1, 3)
                                                   : center;
                    const amrex::Real top = iblank(i, j, k + 1) == -1
                                                ? gpsi(i, j, k + 1, 3)
                                                : center;
                    gpsi(i, j, k, 0) =
                        0.5 * ((right + center) - (center + left)) / dx[0];
                    gpsi(i, j, k, 1) =
                        0.5 * ((front + center) - (center + back)) / dx[1];
                    gpsi(i, j, k, 2) =
                        0.5 * ((top + center) - (center + bottom)) / dx[2];
                } else {
                    gpsi(i, j, k, 0) = 0.;
                    gpsi(i, j, k, 1) = 0.;
                    gpsi(i, j, k, 2) = 0.;
                }
            });
    }
}

void populate_face_vel(
    amrex::MultiFab& mf_vel,
    amrex::iMultiFab& mf_iblank,
    amrex::MultiFab& mf_fx,
    amrex::MultiFab& mf_fy,
    amrex::MultiFab& mf_fz)
{
    // offset for flux indices
    constexpr int ioff = 2;
    // face velocities stored in flux arrays
    for (amrex::MFIter mfi(mf_vel); mfi.isValid(); ++mfi) {
        const auto& vbx = mfi.validbox();
        const auto& xbx = amrex::surroundingNodes(vbx, 0);
        const auto& ybx = amrex::surroundingNodes(vbx, 1);
        // extra points are not really necessary for momentum sharpening
        const auto& zbxp1 = grow(amrex::surroundingNodes(vbx, 2), 1);
        const amrex::Array4<const int>& iblank = mf_iblank.const_array(mfi);
        const amrex::Array4<const amrex::Real>& vel = mf_vel.const_array(mfi);
        const amrex::Array4<amrex::Real>& fx = mf_fx.array(mfi);
        const amrex::Array4<amrex::Real>& fy = mf_fy.array(mfi);
        const amrex::Array4<amrex::Real>& fz = mf_fz.array(mfi);
        // Populate interpolated face velocity
        amrex::ParallelFor(
            xbx, ybx, zbxp1,
            [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                // Neumann condition across nalu bdy
                if (iblank(i, j, k) == -1 || iblank(i - 1, j, k) == -1) {
                    const bool yc = iblank(i, j, k) == -1;
                    const bool yn = iblank(i - 1, j, k) == -1;
                    const amrex::Real ctru = vel(i, j, k, 0);
                    const amrex::Real ctrv = vel(i, j, k, 1);
                    const amrex::Real ctrw = vel(i, j, k, 2);
                    const amrex::Real nbru = vel(i - 1, j, k, 0);
                    const amrex::Real nbrv = vel(i - 1, j, k, 1);
                    const amrex::Real nbrw = vel(i - 1, j, k, 2);
                    fx(i, j, k, ioff + 0) =
                        ((yc ? ctru : 0.0) + (yn ? nbru : 0.0)) /
                        ((yc ? 1.0 : 0.0) + (yn ? 1.0 : 0.0));
                    fx(i, j, k, ioff + 1) =
                        ((yc ? ctrv : 0.0) + (yn ? nbrv : 0.0)) /
                        ((yc ? 1.0 : 0.0) + (yn ? 1.0 : 0.0));
                    fx(i, j, k, ioff + 2) =
                        ((yc ? ctrw : 0.0) + (yn ? nbrw : 0.0)) /
                        ((yc ? 1.0 : 0.0) + (yn ? 1.0 : 0.0));
                } else {
                    fx(i, j, k, ioff + 0) = 0.;
                    fx(i, j, k, ioff + 1) = 0.;
                    fx(i, j, k, ioff + 2) = 0.;
                }
            },
            [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                if (iblank(i, j, k) == -1 || iblank(i, j - 1, k) == -1) {
                    const bool yc = iblank(i, j, k) == -1;
                    const bool yn = iblank(i, j - 1, k) == -1;
                    const amrex::Real ctru = vel(i, j, k, 0);
                    const amrex::Real ctrv = vel(i, j, k, 1);
                    const amrex::Real ctrw = vel(i, j, k, 2);
                    const amrex::Real nbru = vel(i, j - 1, k, 0);
                    const amrex::Real nbrv = vel(i, j - 1, k, 1);
                    const amrex::Real nbrw = vel(i, j - 1, k, 2);
                    fy(i, j, k, ioff + 0) =
                        ((yc ? ctru : 0.0) + (yn ? nbru : 0.0)) /
                        ((yc ? 1.0 : 0.0) + (yn ? 1.0 : 0.0));
                    fy(i, j, k, ioff + 1) =
                        ((yc ? ctrv : 0.0) + (yn ? nbrv : 0.0)) /
                        ((yc ? 1.0 : 0.0) + (yn ? 1.0 : 0.0));
                    fy(i, j, k, ioff + 2) =
                        ((yc ? ctrw : 0.0) + (yn ? nbrw : 0.0)) /
                        ((yc ? 1.0 : 0.0) + (yn ? 1.0 : 0.0));
                } else {
                    fy(i, j, k, ioff + 0) = 0.;
                    fy(i, j, k, ioff + 1) = 0.;
                    fy(i, j, k, ioff + 2) = 0.;
                }
            },
            [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                if (iblank(i, j, k) == -1 || iblank(i, j, k - 1) == -1) {
                    const bool yc = iblank(i, j, k) == -1;
                    const bool yn = iblank(i, j, k - 1) == -1;
                    const amrex::Real ctru = vel(i, j, k, 0);
                    const amrex::Real ctrv = vel(i, j, k, 1);
                    const amrex::Real ctrw = vel(i, j, k, 2);
                    const amrex::Real nbru = vel(i, j, k - 1, 0);
                    const amrex::Real nbrv = vel(i, j, k - 1, 1);
                    const amrex::Real nbrw = vel(i, j, k - 1, 2);
                    fz(i, j, k, ioff + 0) =
                        ((yc ? ctru : 0.0) + (yn ? nbru : 0.0)) /
                        ((yc ? 1.0 : 0.0) + (yn ? 1.0 : 0.0));
                    fz(i, j, k, ioff + 1) =
                        ((yc ? ctrv : 0.0) + (yn ? nbrv : 0.0)) /
                        ((yc ? 1.0 : 0.0) + (yn ? 1.0 : 0.0));
                    fz(i, j, k, ioff + 2) =
                        ((yc ? ctrw : 0.0) + (yn ? nbrw : 0.0)) /
                        ((yc ? 1.0 : 0.0) + (yn ? 1.0 : 0.0));
                } else {
                    fz(i, j, k, ioff + 0) = 0.;
                    fz(i, j, k, ioff + 1) = 0.;
                    fz(i, j, k, ioff + 2) = 0.;
                }
            });
    }
}

amrex::Real alpha_flux(
    int i,
    int j,
    int k,
    int dir,
    const amrex::Real i_th,
    amrex::Array4<amrex::Real const> const& volfrac,
    amrex::Array4<amrex::Real const> const& gdstfnc,
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dx)
{
    // Set up neighbor indices
    int ii = i;
    int jj = j;
    int kk = k;
    ii += (dir == 0) ? -1 : 0;
    jj += (dir == 1) ? -1 : 0;
    kk += (dir == 2) ? -1 : 0;

    // Gradient of phi normal to interface
    const amrex::Real gphi = (volfrac(i, j, k) - volfrac(ii, jj, kk)) / dx[dir];
    // Get gpsi magnitude for average normal
    const amrex::Real gpsi_mag_ = std::sqrt(
        std::pow(gdstfnc(i, j, k, 0), 2) + std::pow(gdstfnc(i, j, k, 1), 2) +
        std::pow(gdstfnc(i, j, k, 2), 2));
    const amrex::Real gpsi_mag_nb = std::sqrt(
        std::pow(gdstfnc(ii, jj, kk, 0), 2) +
        std::pow(gdstfnc(ii, jj, kk, 1), 2) +
        std::pow(gdstfnc(ii, jj, kk, 2), 2));
    const amrex::Real gpsi_norm_avg =
        0.5 * (gdstfnc(ii, jj, kk, dir) / (gpsi_mag_nb + 1e-10) +
               gdstfnc(i, j, k, dir) / (gpsi_mag_ + 1e-10));

    // Determine which psi to use for tanh term
    const amrex::Real psi_upw =
        gphi < 0. ? gdstfnc(i, j, k, 3) : gdstfnc(ii, jj, kk, 3);
    const amrex::Real psi_dnw =
        gphi < 0. ? gdstfnc(ii, jj, kk, 3) : gdstfnc(i, j, k, 3);
    amrex::Real psi_eval = 0.5 * (psi_upw + psi_dnw);
    // Use upwinding to avoid over-emptying
    if ((psi_upw < 0.0 && psi_dnw < 0.0) && psi_dnw > asdf(0.4, i_th)) {
        // Upwind when both upwind and downwind are majority gas, and downwind
        // is not close to interface
        psi_eval = psi_upw;
    }
    const amrex::Real omtanh2_choose =
        1.0 - tanh(0.5 * psi_eval / i_th) * tanh(0.5 * psi_eval / i_th);
    const amrex::Real a_flux =
        (i_th * gphi - 0.25 * omtanh2_choose * gpsi_norm_avg);
    return a_flux;
}

amrex::Real alpha_discrete_flux(
    int i,
    int j,
    int k,
    int dir,
    const amrex::Real margin,
    amrex::Array4<amrex::Real const> const& vof,
    amrex::Array4<amrex::Real const> const& tg_vof,
    amrex::Array4<amrex::Real const> const& normal,
    amrex::Array4<amrex::Real const> const& levelset)
{
    // Set up neighbor indices
    int ii = i;
    int jj = j;
    int kk = k;
    ii += (dir == 0) ? -1 : 0;
    jj += (dir == 1) ? -1 : 0;
    kk += (dir == 2) ? -1 : 0;

    // Gradient of phi normal to interface
    const amrex::Real gphi = (vof(i, j, k) - vof(ii, jj, kk));
    // Find direction of interface using levelset (should be gpsi)
    // const amrex::Real gphi = (levelset(i, j, k) - levelset(ii, jj, kk));
    // Normal vector in each cell (already normalized)
    const amrex::Real norm_ = normal(i, j, k, dir);
    const amrex::Real norm_nb = normal(ii, jj, kk, dir);

    // Determine which delta_phi (and multiply by normal)
    // The sign depends on side of flux face (like upwinding)
    const amrex::Real dphi_ = (tg_vof(i, j, k) - vof(i, j, k)) * (-norm_);
    const amrex::Real dphi_nb =
        (tg_vof(ii, jj, kk) - vof(ii, jj, kk)) * norm_nb;
    // Average value used across the interface
    amrex::Real dphi_eval = 0.5 * (dphi_ + dphi_nb);
    // Upwinding when on the gas side, downwinding on the liquid
    // Across the interface defined as crossing 0.5 or within margin of 0.5
    if ((std::abs(vof(i, j, k) - 0.5) > margin ||
         std::abs(vof(ii, jj, kk) - 0.5) > margin)) {
        if (gphi > 0.0) {
            dphi_eval = (vof(ii, jj, kk) < 0.5 && vof(i, j, k) <= 0.5 + margin)
                            ? dphi_nb
                            : dphi_eval;
            dphi_eval = (vof(ii, jj, kk) >= 0.5 - margin && vof(i, j, k) > 0.5)
                            ? dphi_
                            : dphi_eval;
        }
        if (gphi < 0.0) {
            dphi_eval = (vof(i, j, k) < 0.5 && vof(ii, jj, kk) <= 0.5 + margin)
                            ? dphi_
                            : dphi_eval;
            dphi_eval = (vof(i, j, k) >= 0.5 - margin && vof(ii, jj, kk) > 0.5)
                            ? dphi_nb
                            : dphi_eval;
        }
    }

    /*if (dir == 2 && i == 30 && j == 30 &&
        (std::abs(vof(i, j, k) - 0.5) < 0.1 ||
         std::abs(vof(ii, jj, kk) - 0.5) < 0.1)) {
        std::cout << i << " " << j << " " << k << " " << vof(i, j, k) << " "
                  << vof(ii, jj, kk) << " " << dphi_eval << " " << dphi_ << " "
                  << dphi_nb << std::endl;
    }*/

    return dphi_eval;
}

void velocity_discrete_face(
    int i,
    int j,
    int k,
    int dir,
    const amrex::Real margin,
    amrex::Array4<amrex::Real const> const& vof,
    amrex::Array4<amrex::Real const> const& velocity,
    amrex::Real& uface,
    amrex::Real& vface,
    amrex::Real& wface)
{
    // Set up neighbor indices
    int ii = i;
    int jj = j;
    int kk = k;
    ii += (dir == 0) ? -1 : 0;
    jj += (dir == 1) ? -1 : 0;
    kk += (dir == 2) ? -1 : 0;

    // Gradient of phi normal to interface
    const amrex::Real gphi = (vof(i, j, k) - vof(ii, jj, kk));

    // Get velocities on both sides
    const amrex::Real u_ = velocity(i, j, k, 0);
    const amrex::Real v_ = velocity(i, j, k, 1);
    const amrex::Real w_ = velocity(i, j, k, 2);
    const amrex::Real u_nb = velocity(ii, jj, kk, 0);
    const amrex::Real v_nb = velocity(ii, jj, kk, 1);
    const amrex::Real w_nb = velocity(ii, jj, kk, 2);
    // Average value used across the interface
    uface = 0.5 * (u_ + u_nb);
    vface = 0.5 * (v_ + v_nb);
    wface = 0.5 * (w_ + w_nb);
    // Upwinding when on the gas side, downwinding on the liquid
    // Across the interface defined as crossing 0.5 or within margin of 0.5
    // Same concept as the alpha flux. Get velocity from same cells as density
    if ((std::abs(vof(i, j, k) - 0.5) > margin ||
         std::abs(vof(ii, jj, kk) - 0.5) > margin)) {
        if (gphi > 0.0) {
            bool upw = vof(ii, jj, kk) < 0.5 && vof(i, j, k) <= 0.5 + margin;
            bool dnw = vof(ii, jj, kk) >= 0.5 - margin && vof(i, j, k) > 0.5;
            uface = upw ? u_nb : uface;
            vface = upw ? v_nb : vface;
            wface = upw ? w_nb : wface;
            uface = dnw ? u_ : uface;
            vface = dnw ? v_ : vface;
            wface = dnw ? w_ : wface;
        }
        if (gphi < 0.0) {
            bool upw = vof(i, j, k) < 0.5 && vof(ii, jj, kk) <= 0.5 + margin;
            bool dnw = vof(i, j, k) >= 0.5 - margin && vof(ii, jj, kk) > 0.5;
            uface = upw ? u_ : uface;
            vface = upw ? v_ : vface;
            wface = upw ? w_ : wface;
            uface = dnw ? u_nb : uface;
            vface = dnw ? v_nb : vface;
            wface = dnw ? w_nb : wface;
        }
    }

    /*if (dir == 2 && i == 30 && j == 30 &&
        (std::abs(vof(i, j, k) - 0.5) < 0.1 ||
         std::abs(vof(ii, jj, kk) - 0.5) < 0.1)) {
        std::cout << i << " " << j << " " << k << " " << vof(i, j, k) << " "
                  << vof(ii, jj, kk) << " " << dphi_eval << " " << dphi_ << " "
                  << dphi_nb << std::endl;
    }*/
}

void populate_sharpen_fluxes(
    amrex::MultiFab& mf_fx,
    amrex::MultiFab& mf_fy,
    amrex::MultiFab& mf_fz,
    amrex::MultiFab& mf_gpsi,
    amrex::MultiFab& mf_vof,
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dx,
    const amrex::Real i_th,
    const amrex::Real rho1,
    const amrex::Real rho2)
{
    const amrex::Real dxf = std::cbrt(dx[0] * dx[1] * dx[2]);
    for (amrex::MFIter mfi(mf_vof); mfi.isValid(); ++mfi) {
        const auto& vbx = mfi.validbox();
        const auto& xbx = amrex::surroundingNodes(vbx, 0);
        const auto& ybx = amrex::surroundingNodes(vbx, 1);
        // Extra points for zbox are for pressure source
        const auto& zbxp1 = grow(amrex::surroundingNodes(vbx, 2), 1);
        const amrex::Array4<amrex::Real>& fx = mf_fx.array(mfi);
        const amrex::Array4<amrex::Real>& fy = mf_fy.array(mfi);
        const amrex::Array4<amrex::Real>& fz = mf_fz.array(mfi);
        const amrex::Array4<const amrex::Real>& gpsi = mf_gpsi.const_array(mfi);
        const amrex::Array4<const amrex::Real>& vof = mf_vof.const_array(mfi);
        // Populate vof and density fluxes for each direction
        amrex::ParallelFor(
            xbx, ybx, zbxp1,
            [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                // vof flux
                amrex::Real flux = alpha_flux(i, j, k, 0, i_th, vof, gpsi, dx);
                fx(i, j, k, 0) = flux;
                // density flux
                flux *= (rho1 - rho2);
                fx(i, j, k, 1) = flux;
                // momentum fluxes (dens flux * face vel)
                fx(i, j, k, 2) *= flux;
                fx(i, j, k, 3) *= flux;
                fx(i, j, k, 4) *= flux;
            },
            [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                amrex::Real flux = alpha_flux(i, j, k, 1, i_th, vof, gpsi, dx);
                fy(i, j, k, 0) = flux;
                flux *= (rho1 - rho2);
                fy(i, j, k, 1) = flux;
                fy(i, j, k, 2) *= flux;
                fy(i, j, k, 3) *= flux;
                fy(i, j, k, 4) *= flux;
            },
            [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                amrex::Real flux = alpha_flux(i, j, k, 2, i_th, vof, gpsi, dx);
                fz(i, j, k, 0) = flux;
                flux *= (rho1 - rho2);
                fz(i, j, k, 1) = flux;
                fz(i, j, k, 2) *= flux;
                fz(i, j, k, 3) *= flux;
                fz(i, j, k, 4) *= flux;
                // Turn "on" all z faces, later modified in process_fluxes
                fz(i, j, k, 5) = 1.0;
            });
    }
}

void populate_sharpen_discrete_fluxes(
    amrex::MultiFab& mf_fx,
    amrex::MultiFab& mf_fy,
    amrex::MultiFab& mf_fz,
    amrex::MultiFab& mf_vof,
    amrex::MultiFab& mf_target_vof,
    amrex::MultiFab& mf_norm,
    amrex::MultiFab& mf_levelset,
    amrex::MultiFab& mf_velocity,
    const amrex::Real margin,
    const amrex::Real rho1,
    const amrex::Real rho2)
{
    for (amrex::MFIter mfi(mf_vof); mfi.isValid(); ++mfi) {
        const auto& vbx = mfi.validbox();
        const auto& xbx = amrex::surroundingNodes(vbx, 0);
        const auto& ybx = amrex::surroundingNodes(vbx, 1);
        // Extra points for zbox are for pressure source
        const auto& zbxp1 = grow(amrex::surroundingNodes(vbx, 2), 1);
        const amrex::Array4<amrex::Real>& fx = mf_fx.array(mfi);
        const amrex::Array4<amrex::Real>& fy = mf_fy.array(mfi);
        const amrex::Array4<amrex::Real>& fz = mf_fz.array(mfi);
        const amrex::Array4<const amrex::Real>& vof = mf_vof.const_array(mfi);
        const amrex::Array4<const amrex::Real>& tg_vof =
            mf_target_vof.const_array(mfi);
        const amrex::Array4<const amrex::Real>& norm = mf_norm.const_array(mfi);
        const amrex::Array4<const amrex::Real>& lvs =
            mf_levelset.const_array(mfi);
        const amrex::Array4<const amrex::Real>& vel =
            mf_velocity.const_array(mfi);
        // Populate vof and density fluxes for each direction
        amrex::ParallelFor(
            xbx, ybx, zbxp1,
            [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                // vof flux
                amrex::Real flux = alpha_discrete_flux(
                    i, j, k, 0, margin, vof, tg_vof, norm, lvs);
                fx(i, j, k, 0) = flux;
                // density flux
                flux *= (rho1 - rho2);
                fx(i, j, k, 1) = flux;
                // momentum fluxes (dens flux * face vel)
                amrex::Real uf, vf, wf;
                velocity_discrete_face(
                    i, j, k, 0, margin, vof, vel, uf, vf, wf);
                fx(i, j, k, 2) = flux * uf;
                fx(i, j, k, 3) = flux * vf;
                fx(i, j, k, 4) = flux * wf;
            },
            [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                amrex::Real flux = alpha_discrete_flux(
                    i, j, k, 1, margin, vof, tg_vof, norm, lvs);
                fy(i, j, k, 0) = flux;
                flux *= (rho1 - rho2);
                fy(i, j, k, 1) = flux;
                amrex::Real uf, vf, wf;
                velocity_discrete_face(
                    i, j, k, 1, margin, vof, vel, uf, vf, wf);
                fy(i, j, k, 2) = flux * uf;
                fy(i, j, k, 3) = flux * vf;
                fy(i, j, k, 4) = flux * wf;
            },
            [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                amrex::Real flux = alpha_discrete_flux(
                    i, j, k, 2, margin, vof, tg_vof, norm, lvs);
                fz(i, j, k, 0) = flux;
                flux *= (rho1 - rho2);
                fz(i, j, k, 1) = flux;
                amrex::Real uf, vf, wf;
                velocity_discrete_face(
                    i, j, k, 2, margin, vof, vel, uf, vf, wf);
                fz(i, j, k, 2) = flux * uf;
                fz(i, j, k, 3) = flux * vf;
                fz(i, j, k, 4) = flux * wf;
                // Turn "on" all z faces, later modified in process_fluxes
                fz(i, j, k, 5) = 1.0;

                /*if (i == 30 && j == 30 && (k == 31 || k == 32 || k == 33)) {
                    std::cout << "populate_fluxes " << i << " " << j << " " << k
                              << " " << fz(i, j, k, 0) << std::endl;
                }*/
            });
    }
}

void process_fluxes_calc_src(
    amrex::MultiFab& mf_fx,
    amrex::MultiFab& mf_fy,
    amrex::MultiFab& mf_fz,
    amrex::MultiFab& mf_psource,
    amrex::iMultiFab& mf_iblank,
    const amrex::Real grav_z)
{
    int ncompx = mf_fx.nComp();
    int ncompy = mf_fy.nComp();
    int ncompz = mf_fz.nComp();
    for (amrex::MFIter mfi(mf_iblank); mfi.isValid(); ++mfi) {
        const auto& vbx = mfi.validbox();
        const auto& xbx = amrex::surroundingNodes(vbx, 0);
        const auto& ybx = amrex::surroundingNodes(vbx, 1);
        const auto& zbxp1 = grow(amrex::surroundingNodes(vbx, 2), 1);
        const amrex::Array4<amrex::Real>& fx = mf_fx.array(mfi);
        const amrex::Array4<amrex::Real>& fy = mf_fy.array(mfi);
        const amrex::Array4<amrex::Real>& fz = mf_fz.array(mfi);
        const amrex::Array4<amrex::Real>& sp = mf_psource.array(mfi);
        const amrex::Array4<const int>& iblank = mf_iblank.const_array(mfi);
        // Only faces with iblank = -1 on both sides can have nonzero flux
        amrex::ParallelFor(
            xbx, ybx, zbxp1,
            [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                bool zero_all = (iblank(i - 1, j, k) + iblank(i, j, k) > -2);
                for (int n = 0; n < ncompx; ++n) {
                    fx(i, j, k, n) *= zero_all ? 0. : 1.;
                }
            },
            [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                bool zero_all = (iblank(i, j - 1, k) + iblank(i, j, k) > -2);
                for (int n = 0; n < ncompy; ++n) {
                    fy(i, j, k, n) *= zero_all ? 0. : 1.;
                }
            },
            [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                bool zero_all = (iblank(i, j, k - 1) + iblank(i, j, k) > -2);
                for (int n = 0; n < ncompz; ++n) {
                    fz(i, j, k, n) *= zero_all ? 0. : 1.;
                }
                /*if (i == 30 && j == 30 && (k == 31 || k == 32 || k == 33)) {
                    std::cout
                        << "process_fluxes " << i << " " << j << " " << k << " "
                        << fz(i, j, k, 0) << " zero_all " << zero_all << " "
                        << iblank(i, j, k - 1) + iblank(i, j, k) << std::endl;
                }*/
            });
        // With knowledge of fluxes, compute pressure source term
        amrex::Box const& nbx = mfi.nodaltilebox();
        amrex::ParallelFor(
            nbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                // To get the local density flux corresponding to the node,
                // average the surrounding fluxes. Fluxes in AMR-Wind domain are
                // zeroed, and they are flagged as 0 to avoid dividing by the
                // wrong number in the average. For AMR-Wind nodes (not
                // bordering nalu-wind cells), psource = 0 / (0 + tiny)
                const amrex::Real rho_flux =
                    (fz(i, j, k, 1) + fz(i - 1, j, k, 1) + fz(i, j - 1, k, 1) +
                     fz(i - 1, j - 1, k, 1)) /
                    (fz(i, j, k, 5) + fz(i - 1, j, k, 5) + fz(i, j - 1, k, 5) +
                     fz(i - 1, j - 1, k, 5) + 1e-100);
                // Pressure source is proportional to local density source
                // (because more or less density is above node)
                sp(i, j, k) = rho_flux * grav_z;
            });
    }
}

void apply_fluxes(
    amrex::MultiFab& mf_fx,
    amrex::MultiFab& mf_fy,
    amrex::MultiFab& mf_fz,
    amrex::MultiFab& mf_psource,
    amrex::MultiFab& mf_vof,
    amrex::MultiFab& mf_dens,
    amrex::MultiFab& mf_vel,
    amrex::MultiFab& mf_pressure,
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dx,
    const amrex::Real i_th)
{
    int ncomp = mf_fx.nComp();
    constexpr amrex::Real tiny = 1e-12;
    // Pseudo-time factor
    const amrex::Real ptfac =
        i_th > 0. ? 1.5 * i_th : std::cbrt(dx[0] * dx[1] * dx[2]);
    for (amrex::MFIter mfi(mf_vof); mfi.isValid(); ++mfi) {
        const auto& vbx = mfi.validbox();
        amrex::Box const& nbx = mfi.nodaltilebox();
        const amrex::Array4<amrex::Real>& fx = mf_fx.array(mfi);
        const amrex::Array4<amrex::Real>& fy = mf_fy.array(mfi);
        const amrex::Array4<amrex::Real>& fz = mf_fz.array(mfi);
        const amrex::Array4<amrex::Real>& vof = mf_vof.array(mfi);
        const amrex::Array4<amrex::Real>& dens = mf_dens.array(mfi);
        const amrex::Array4<amrex::Real>& vel = mf_vel.array(mfi);
        const amrex::Array4<amrex::Real>& p = mf_pressure.array(mfi);
        const amrex::Array4<amrex::Real>& sp = mf_psource.array(mfi);
        amrex::ParallelFor(
            vbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                const amrex::Real oldvof = vof(i, j, k);
                const amrex::Real olddens = dens(i, j, k);
                vof(i, j, k) +=
                    ptfac * ((fx(i + 1, j, k, 0) - fx(i, j, k, 0)) / dx[0] +
                             (fy(i, j + 1, k, 0) - fy(i, j, k, 0)) / dx[1] +
                             (fz(i, j, k + 1, 0) - fz(i, j, k, 0)) / dx[2]);
                dens(i, j, k) +=
                    ptfac * ((fx(i + 1, j, k, 1) - fx(i, j, k, 1)) / dx[0] +
                             (fy(i, j + 1, k, 1) - fy(i, j, k, 1)) / dx[1] +
                             (fz(i, j, k + 1, 1) - fz(i, j, k, 1)) / dx[2]);
                vel(i, j, k, 0) =
                    1.0 / dens(i, j, k) *
                    (olddens * vel(i, j, k, 0) +
                     ptfac * ((fx(i + 1, j, k, 2) - fx(i, j, k, 2)) / dx[0] +
                              (fy(i, j + 1, k, 2) - fy(i, j, k, 2)) / dx[1] +
                              (fz(i, j, k + 1, 2) - fz(i, j, k, 2)) / dx[2]));
                vel(i, j, k, 1) =
                    1.0 / dens(i, j, k) *
                    (olddens * vel(i, j, k, 1) +
                     ptfac * ((fx(i + 1, j, k, 3) - fx(i, j, k, 3)) / dx[0] +
                              (fy(i, j + 1, k, 3) - fy(i, j, k, 3)) / dx[1] +
                              (fz(i, j, k + 1, 3) - fz(i, j, k, 3)) / dx[2]));
                vel(i, j, k, 2) =
                    1.0 / dens(i, j, k) *
                    (olddens * vel(i, j, k, 2) +
                     ptfac * ((fx(i + 1, j, k, 4) - fx(i, j, k, 4)) / dx[0] +
                              (fy(i, j + 1, k, 4) - fy(i, j, k, 4)) / dx[1] +
                              (fz(i, j, k + 1, 4) - fz(i, j, k, 4)) / dx[2]));

                /*if (i == 30 && j == 30 && (k == 31 || k == 32)) {
                    std::cout << "apply_fluxes " << i << " " << j << " " << k
                              << " " << oldvof << " " << vof(i, j, k) << " "
                              << fz(i, j, k + 1, 0) << " " << -fz(i, j, k, 0)
                              << std::endl;
                }*/

                // Ensure vof is bounded
                vof(i, j, k) =
                    vof(i, j, k) < tiny
                        ? 0.0
                        : (vof(i, j, k) > 1. - tiny ? 1. : vof(i, j, k));
                // Density is corrected later
            });

        // Apply pressure source/sink
        amrex::ParallelFor(
            nbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                // Move hydrostatic pressure around due to density flux
                p(i, j, k) += ptfac * sp(i, j, k);
            });
    }
}

amrex::Real measure_convergence(
    amrex::MultiFab& mf_fx, amrex::MultiFab& mf_fy, amrex::MultiFab& mf_fz)
{
    // Get the maximum flux magnitude, but just for vof fluxes
    const amrex::Real err_fx = amrex::ReduceMax(
        mf_fx, 0,
        [=] AMREX_GPU_HOST_DEVICE(
            amrex::Box const& bx,
            amrex::Array4<amrex::Real const> const& fx) -> amrex::Real {
            amrex::Real err_fab = -1.0;
            amrex::Loop(bx, [=, &err_fab](int i, int j, int k) noexcept {
                err_fab = amrex::max(err_fab, std::abs(fx(i, j, k, 0)));
            });
            return err_fab;
        });
    const amrex::Real err_fy = amrex::ReduceMax(
        mf_fy, 0,
        [=] AMREX_GPU_HOST_DEVICE(
            amrex::Box const& bx,
            amrex::Array4<amrex::Real const> const& fy) -> amrex::Real {
            amrex::Real err_fab = -1.0;
            amrex::Loop(bx, [=, &err_fab](int i, int j, int k) noexcept {
                err_fab = amrex::max(err_fab, std::abs(fy(i, j, k, 0)));
            });
            return err_fab;
        });
    const amrex::Real err_fz = amrex::ReduceMax(
        mf_fz, 0,
        [=] AMREX_GPU_HOST_DEVICE(
            amrex::Box const& bx,
            amrex::Array4<amrex::Real const> const& fz) -> amrex::Real {
            amrex::Real err_fab = -1.0;
            amrex::Loop(bx, [=, &err_fab](int i, int j, int k) noexcept {
                err_fab = amrex::max(err_fab, std::abs(fz(i, j, k, 0)));
            });
            return err_fab;
        });
    const amrex::Real err = amrex::max(err_fx, amrex::max(err_fy, err_fz));
    return err;
}

// for debugging
void equate_vof(amrex::MultiFab& mf_gpsi, amrex::MultiFab& mf_vof)
{
    for (amrex::MFIter mfi(mf_vof); mfi.isValid(); ++mfi) {
        const auto& vbx = mfi.validbox();
        const amrex::Array4<const amrex::Real>& gpsi = mf_gpsi.const_array(mfi);
        const amrex::Array4<amrex::Real>& vof = mf_vof.array(mfi);
        amrex::ParallelFor(
            vbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                vof(i, j, k) = gpsi(i, j, k, 0);
            });
    }
}

void harmonize_vof(
    amrex::MultiFab& mf_vof_target,
    amrex::MultiFab& mf_vof_original,
    amrex::iMultiFab& mf_iblank)
{
    for (amrex::MFIter mfi(mf_vof_target); mfi.isValid(); ++mfi) {
        const auto& vbx = mfi.validbox();
        const amrex::Array4<amrex::Real>& tg_vof = mf_vof_target.array(mfi);
        const amrex::Array4<const amrex::Real>& og_vof =
            mf_vof_original.const_array(mfi);
        const amrex::Array4<const int>& iblank = mf_iblank.const_array(mfi);
        amrex::ParallelFor(
            vbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                // Replace amr-wind vof values with originals
                if (iblank(i, j, k) != -1) {
                    tg_vof(i, j, k) = og_vof(i, j, k);
                }
            });
    }
}

void SharpenNaluData(
    CFDSim& sim,
    const int nit,
    const amrex::Real tol,
    const int ncconv,
    const amrex::Real rlscale)
{
    auto& repo = sim.repo();
    auto nlevels = repo.num_active_levels();
    auto geom = sim.mesh().Geom();

    // Check for vof variable -- sharpening only required with VOF
    bool yesVOF = repo.field_exists("vof");
    if (!yesVOF) return;

    // Get phase densities
    auto& mphase = sim.physics_manager().get<MultiPhase>();
    const amrex::Real rhol = mphase.rho1();
    const amrex::Real rhog = mphase.rho2();
    const amrex::Real grav_z = mphase.gravity()[2];

    // Get blanking for cells
    auto& iblank_cell = repo.get_int_field("iblank_cell");

    // Get fields that will be modified
    auto& vof = repo.get_field("vof");
    auto& rho = repo.get_field("density");
    auto& velocity = repo.get_field("velocity");
    auto& p = repo.get_field("p");

    // Create scratch fields for fluxes - 5 components are vof, density, and
    // 3 of velocity (flux_z has 6th component for flag in p_src calc)
    auto flux_x = repo.create_scratch_field(5, 0, amr_wind::FieldLoc::XFACE);
    auto flux_y = repo.create_scratch_field(5, 0, amr_wind::FieldLoc::YFACE);
    auto flux_z = repo.create_scratch_field(6, 1, amr_wind::FieldLoc::ZFACE);
    // Create scratch field for pressure source term
    auto p_src = repo.create_scratch_field(1, 0, amr_wind::FieldLoc::NODE);
    // Create scratch field for approximate signed distance function and grad
    // (components 0-2 are gradient, 3 is asdf)
    auto grad_asdf = repo.create_scratch_field(4, vof.num_grow()[0]);

    // Process vof for tiny margins from single-phase
    /*for (int lev = 0; lev < nlevels; ++lev) {
        process_vof(vof(lev));
    }*/

    // Pseudo-time loop
    amrex::Real err = 100.0 * tol;
    int n = 0;
    while (n < nit && err > tol) {
        // Increment step counter
        ++n;
        // Determine if convergence error is calculated this step
        bool cconv = n % ncconv == 0;
        // Zero error if being calculated this step
        err = cconv ? 0.0 : err;

        for (int lev = 0; lev < nlevels; ++lev) {
            // Thickness of interface should be 1/2 the size of the cell?
            auto dx = (geom[lev]).CellSizeArray();
            const amrex::Real i_th = rlscale * std::cbrt(dx[0] * dx[1] * dx[2]);

            // Populate approximate signed distance function
            populate_grad_psi(
                (*grad_asdf)(lev), vof(lev), iblank_cell(lev), dx, i_th);

            // Populate interpolated face velocities, store in flux arrays
            populate_face_vel(
                velocity(lev), iblank_cell(lev), (*flux_x)(lev), (*flux_y)(lev),
                (*flux_z)(lev));

            // Sharpening fluxes for vof, density, pressure gradient
            populate_sharpen_fluxes(
                (*flux_x)(lev), (*flux_y)(lev), (*flux_z)(lev),
                (*grad_asdf)(lev), vof(lev), dx, i_th, rhol, rhog);

            // Zero fluxes outside nalu portion of domain and calculate p_src
            process_fluxes_calc_src(
                (*flux_x)(lev), (*flux_y)(lev), (*flux_z)(lev), (*p_src)(lev),
                iblank_cell(lev), grav_z);

            // Apply fluxes
            apply_fluxes(
                (*flux_x)(lev), (*flux_y)(lev), (*flux_z)(lev), (*p_src)(lev),
                vof(lev), rho(lev), velocity(lev), p(lev), dx, i_th);

            // Measure convergence to determine if loop can stop
            if (cconv) {
                // Update error at specified interval of steps
                const amrex::Real err_lev = measure_convergence(
                    (*flux_x)(lev), (*flux_y)(lev), (*flux_z)(lev));
                err = amrex::max(err, err_lev);
            }

            // Replace vof with another field to look at it
            // equate_vof((*grad_asdf)(lev), vof(lev));
        }

        // Average down fluxes
        // Get pseudo CFL (multiplier for sharpening term)
        // Apply fluxes

        // Average down vof for consistency

        // Fillpatch for ghost cells
        vof.fillpatch(sim.time().current_time());
        velocity.fillpatch(sim.time().current_time());

        // Update density (fillpatch built in)
        mphase.set_density_via_vof();

        // Ensure that err is same across processors
        if (cconv) {
            amrex::ParallelDescriptor::ReduceRealMax(err);
        }

        amrex::Print() << "sharpen step " << n << " " << err << " " << tol
                       << std::endl;
    }

    // Copy to old
    amr_wind::field_ops::copy(
        vof.state(amr_wind::FieldState::Old), vof, 0, 0, vof.num_comp(),
        vof.num_grow());
    amr_wind::field_ops::copy(
        rho.state(amr_wind::FieldState::Old), rho, 0, 0, rho.num_comp(),
        rho.num_grow());
    amr_wind::field_ops::copy(
        velocity.state(amr_wind::FieldState::Old), velocity, 0, 0,
        velocity.num_comp(), velocity.num_grow());

    // Fillpatch for pressure to make sure pressure stencil has all points
    p.fillpatch(sim.time().current_time());
}

void SharpenNaluDataDiscrete(
    CFDSim& sim,
    const int nit,
    const amrex::Real tol,
    const int ncconv,
    const amrex::Real rlscale,
    const amrex::Real margin,
    const bool sh_pressure)
{
    auto& repo = sim.repo();
    auto nlevels = repo.num_active_levels();
    auto geom = sim.mesh().Geom();

    // Check for vof variable -- sharpening only required with VOF
    bool yesVOF = repo.field_exists("vof");
    if (!yesVOF) return;

    // Get phase densities
    auto& mphase = sim.physics_manager().get<MultiPhase>();
    const amrex::Real rhol = mphase.rho1();
    const amrex::Real rhog = mphase.rho2();
    // If pressure is not to be modified, set gravity to 0
    const amrex::Real grav_z = sh_pressure ? mphase.gravity()[2] : 0.0;

    // Get blanking for cells
    auto& iblank_cell = repo.get_int_field("iblank_cell");

    // Get fields that will be modified
    auto& vof = repo.get_field("vof");
    auto& levelset = repo.get_field("levelset");
    auto& rho = repo.get_field("density");
    auto& velocity = repo.get_field("velocity");
    auto& p = repo.get_field("p");

    // Create scratch fields for fluxes - 5 components are vof, density, and
    // 3 of velocity (flux_z has 6th component for flag in p_src calc)
    auto flux_x = repo.create_scratch_field(5, 0, amr_wind::FieldLoc::XFACE);
    auto flux_y = repo.create_scratch_field(5, 0, amr_wind::FieldLoc::YFACE);
    auto flux_z = repo.create_scratch_field(6, 1, amr_wind::FieldLoc::ZFACE);
    // Create scratch field for pressure source term
    auto p_src = repo.create_scratch_field(1, 0, amr_wind::FieldLoc::NODE);
    // Create scratch field for approximate signed distance function and grad
    // (components 0-2 are gradient, 3 is asdf)
    auto normal_vec = repo.create_scratch_field(3, vof.num_grow()[0] - 1);

    auto target_vof = repo.create_scratch_field(1, vof.num_grow()[0]);

    // Process vof for tiny margins from single-phase
    /*for (int lev = 0; lev < nlevels; ++lev) {
        process_vof(vof(lev));
    }*/

    // Create levelset field
    for (int lev = 0; lev < nlevels; ++lev) {
        // Thickness used here is user parameter, whatever works best
        auto dx = (geom[lev]).CellSizeArray();
        const amrex::Real i_th = rlscale * std::cbrt(dx[0] * dx[1] * dx[2]);

        // Populate approximate signed distance function
        populate_psi(levelset(lev), vof(lev), i_th);
    }

    // Convert levelset to vof to get target_vof
    mphase.levelset2vof(iblank_cell, *target_vof);

    // Replace vof with original values in amr domain
    for (int lev = 0; lev < nlevels; ++lev) {
        harmonize_vof((*target_vof)(lev), vof(lev), iblank_cell(lev));
    }

    // Pseudo-time loop
    amrex::Real err = 100.0 * tol;
    int n = 0;
    while (n < nit && err > tol) {
        // Increment step counter
        ++n;
        // Determine if convergence error is calculated this step
        bool cconv = n % ncconv == 0;
        // Zero error if being calculated this step
        err = cconv ? 0.0 : err;

        for (int lev = 0; lev < nlevels; ++lev) {
            // Thickness used here is user parameter, whatever works best
            auto dx = (geom[lev]).CellSizeArray();
            // const amrex::Real i_th = rlscale * std::cbrt(dx[0] * dx[1] *
            // dx[2]);

            // Populate normal vector
            populate_normal_vector(
                (*normal_vec)(lev), vof(lev), iblank_cell(lev));

            // Populate interpolated face velocities, store in flux arrays
            /*populate_face_vel(
                velocity(lev), iblank_cell(lev), (*flux_x)(lev), (*flux_y)(lev),
                (*flux_z)(lev));*/

            // Sharpening fluxes for vof, density, and momentum
            populate_sharpen_discrete_fluxes(
                (*flux_x)(lev), (*flux_y)(lev), (*flux_z)(lev), vof(lev),
                (*target_vof)(lev), (*normal_vec)(lev), levelset(lev),
                velocity(lev), margin, rhol, rhog);

            // Process fluxes and get pressure source term
            process_fluxes_calc_src(
                (*flux_x)(lev), (*flux_y)(lev), (*flux_z)(lev), (*p_src)(lev),
                iblank_cell(lev), grav_z);

            // Apply fluxes
            apply_fluxes(
                (*flux_x)(lev), (*flux_y)(lev), (*flux_z)(lev), (*p_src)(lev),
                vof(lev), rho(lev), velocity(lev), p(lev), dx, -1.0);

            // Measure convergence to determine if loop can stop
            if (cconv) {
                // Update error at specified interval of steps
                const amrex::Real err_lev = measure_convergence(
                    (*flux_x)(lev), (*flux_y)(lev), (*flux_z)(lev));
                err = amrex::max(err, err_lev);
            }
        }

        // sim.io_manager().write_plot_file();
        // sim.time().increment_timestep();

        // Average down fluxes
        // Get pseudo CFL (multiplier for sharpening term)
        // Apply fluxes

        // Average down vof for consistency

        // Fillpatch for ghost cells
        vof.fillpatch(sim.time().current_time());
        velocity.fillpatch(sim.time().current_time());

        // Update density (fillpatch built in)
        mphase.set_density_via_vof();

        // Ensure that err is same across processors
        if (cconv) {
            amrex::ParallelDescriptor::ReduceRealMax(err);
        }

        amrex::Print() << "sharpen step " << n << " " << err << " " << tol
                       << std::endl;

        // std::cout << n << " new step\n";
    }

    /*for (int lev = 0; lev < nlevels; ++lev) {
        equate_vof((*target_vof)(lev), vof(lev));
    }*/

    // Copy to old
    amr_wind::field_ops::copy(
        vof.state(amr_wind::FieldState::Old), vof, 0, 0, vof.num_comp(),
        vof.num_grow());
    amr_wind::field_ops::copy(
        rho.state(amr_wind::FieldState::Old), rho, 0, 0, rho.num_comp(),
        rho.num_grow());
    amr_wind::field_ops::copy(
        velocity.state(amr_wind::FieldState::Old), velocity, 0, 0,
        velocity.num_comp(), velocity.num_grow());

    // Fillpatch for pressure to make sure pressure stencil has all points
    p.fillpatch(sim.time().current_time());
}
} // namespace overset
} // namespace amr_wind

#endif