#ifndef SHARPEN_NALU_DATA_H_
#define SHARPEN_NALU_DATA_H_

#include "amr-wind/CFDSim.H"
#include "amr-wind/physics/multiphase/MultiPhase.H"
#include "amr-wind/equation_systems/vof/volume_fractions.H"

namespace amr_wind {
namespace overset {

// Approximate signed distance function
amrex::Real asdf(const amrex::Real a_vof, const amrex::Real i_th)
{
    // function of local vof value and interface thickness
    return (i_th * log((a_vof + 1e-100) / (1. - a_vof + 1e-100)));
}

// Hyperbolic tangent (squared) with a limit
amrex::Real tanh2_lim(const amrex::Real arg)
{
    constexpr amrex::Real tol = 1e-12;
    amrex::Real result = tanh(arg) * tanh(arg);
    // Force tanh result based on limits of asdf calc
    if (arg < 0.5 * log(1e-100 / (1.0 + 1e-100)) + tol) {
        result = 1.0;
    } else if (arg > 0.5 * log((1.0 + 1e-100) / 1e-100) - tol) {
        result = 1.0;
    }
    return result;
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

amrex::Real alpha_flux(
    int i,
    int j,
    int k,
    int dir,
    const amrex::Real i_th,
    amrex::Array4<amrex::Real const> const& volfrac,
    amrex::Array4<amrex::Real const> const& gdstfnc,
    amrex::Array4<int const> const& iblank,
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dx)
{
    // Set up neighbor indices
    int ii = i;
    int jj = j;
    int kk = k;
    ii += (dir == 0) ? -1 : 0;
    jj += (dir == 1) ? -1 : 0;
    kk += (dir == 2) ? -1 : 0;

    // Get face quantities
    // Gradient of phi normal to interface
    const amrex::Real gphi = (volfrac(i, j, k) - volfrac(ii, jj, kk)) / dx[dir];
    const amrex::Real psi_avg =
        0.5 * (gdstfnc(i, j, k, 3) + gdstfnc(ii, jj, kk, 3));
    // Get average gpsi magnitude, replace
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
    // Determine bcs for each cell
    int ibdy_ = 0, jbdy_ = 0, kbdy_ = 0;
    int ibdy_nb = 0, jbdy_nb = 0, kbdy_nb = 0;
    ibdy_ += (iblank(i, j, k) + iblank(i - 1, j, k) > -2)
                 ? -1
                 : (iblank(i, j, k) + iblank(i + 1, j, k) > -2 ? +1 : 0);
    ibdy_nb +=
        (iblank(ii, jj, k) + iblank(ii - 1, jj, kk) > -2)
            ? -1
            : (iblank(ii, jj, kk) + iblank(ii + 1, jj, kk) > -2 ? +1 : 0);
    jbdy_ += (iblank(i, j, k) + iblank(i, j - 1, k) > -2)
                 ? -1
                 : (iblank(i, j, k) + iblank(i, j + 1, k) > -2 ? +1 : 0);
    jbdy_nb +=
        (iblank(ii, jj, k) + iblank(ii, jj - 1, kk) > -2)
            ? -1
            : (iblank(ii, jj, kk) + iblank(ii, jj + 1, kk) > -2 ? +1 : 0);
    kbdy_ += (iblank(i, j, k) + iblank(i, j, k - 1) > -2)
                 ? -1
                 : (iblank(i, j, k) + iblank(i, j, k + 1) > -2 ? +1 : 0);
    kbdy_nb +=
        (iblank(ii, jj, k) + iblank(ii, jj, kk - 1) > -2)
            ? -1
            : (iblank(ii, jj, kk) + iblank(ii, jj, kk + 1) > -2 ? +1 : 0);
    // Get normals for each cell
    amrex::Real mx_ = 0., my_ = 0., mz_ = 0.;
    amrex::Real mx_nb = 0., my_nb = 0., mz_nb = 0.;
    amr_wind::multiphase::youngs_fd_normal_neumann(
        i, j, k, ibdy_, jbdy_, kbdy_, volfrac, mx_, my_, mz_);
    amr_wind::multiphase::youngs_fd_normal_neumann(
        ii, jj, kk, ibdy_nb, jbdy_nb, kbdy_nb, volfrac, mx_nb, my_nb, mz_nb);
    amrex::Real mn_ = dir == 0 ? mx_ : (dir == 1 ? my_ : mz_);
    const amrex::Real mn_nb = dir == 0 ? mx_nb : (dir == 1 ? my_nb : mz_nb);
    // Average normals from each and divide by dx (not done by other)
    mx_ = 0.5 * (mx_ + mx_nb) / dx[0];
    my_ = 0.5 * (my_ + my_nb) / dx[1];
    mz_ = 0.5 * (mz_ + mz_nb) / dx[2];
    mn_ = 0.5 * (mn_ + mn_nb) / dx[dir];
    // Calculate magnitude and normalize, replacing normal gradient
    const amrex::Real gphi_norm = gphi / std::sqrt(
                                             mx_ * mx_ + my_ * my_ + mz_ * mz_ -
                                             mn_ * mn_ + gphi * gphi + 1e-100);

    // Tanh term - using average of tanh2 in each cell
    amrex::Real omtanh2_avg =
        1.0 - 0.5 * (tanh2_lim(0.5 * gdstfnc(i, j, k, 3) / i_th) +
                     tanh2_lim(0.5 * gdstfnc(ii, jj, kk, 3) / i_th));
    amrex::Real omtanh_avg2 =
        1.0 - tanh(0.5 * psi_avg / i_th) * tanh(0.5 * psi_avg / i_th);
    const amrex::Real psi_upw =
        gphi < 0. ? gdstfnc(i, j, k, 3) : gdstfnc(ii, jj, kk, 3);
    const amrex::Real psi_dnw =
        gphi < 0. ? gdstfnc(ii, jj, kk, 3) : gdstfnc(i, j, k, 3);
    amrex::Real psi_eval = 0.5 * (psi_upw + psi_dnw);
    if (psi_upw * psi_dnw >= 0.0) {
        // Upwind and downwind are same majority
        if (psi_upw < 0.0) {
            // Upwind in majority gas
            psi_eval = psi_upw;
        }
        if (psi_dnw > 0.0) {
            // Downwind in majority liquid
            psi_eval = psi_dnw;
        }
    }
    const amrex::Real omtanh2_upw = 1.0 - tanh2_lim(0.5 * psi_upw / i_th);
    const amrex::Real omtanh2_choose = 1.0 - tanh2_lim(0.5 * psi_eval / i_th);

    /*if (i == 20 && j == 20 && omtanh_avg2 > 0.0 &&
        std::abs(gphi_norm) > 1e-12) {
        std::cout << "omtanh " << i << " " << j << " " << k << " " << dir << " "
                  << omtanh_avg2 << " " << omtanh2_avg << std::endl;
    }*/

    const amrex::Real a_flux =
        (i_th * gphi - 0.25 * omtanh2_choose * gphi_norm);
    // const amrex::Real a_flux = (i_th * gphi - 0.25 * omtanh2_upw *
    // gphi_norm);
    //  const amrex::Real a_flux = (i_th * gphi - 0.25 * omtanh2_avg *
    //  gphi_norm); const amrex::Real a_flux = (i_th * gphi - 0.25 * omtanh_avg2
    //  * gphi_norm); const amrex::Real a_flux = (i_th * gphi - 0.25 *
    //  omtanh_avg2
    //  * gpsi_norm_avg); const amrex::Real a_flux = (- 0.25 * omtanh_avg2 *
    //  gpsi_norm_avg);

    // Would be multiplied by Gamma (velocity scale)
    return a_flux;
}

void populate_sharpen_fluxes(
    amrex::MultiFab& mf_fx,
    amrex::MultiFab& mf_fy,
    amrex::MultiFab& mf_fz,
    amrex::MultiFab& mf_gpsi,
    amrex::MultiFab& mf_vof,
    amrex::iMultiFab& mf_iblank,
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
        const auto& zbx = amrex::surroundingNodes(vbx, 2);
        const amrex::Array4<amrex::Real>& fx = mf_fx.array(mfi);
        const amrex::Array4<amrex::Real>& fy = mf_fy.array(mfi);
        const amrex::Array4<amrex::Real>& fz = mf_fz.array(mfi);
        const amrex::Array4<const amrex::Real>& gpsi = mf_gpsi.const_array(mfi);
        const amrex::Array4<const amrex::Real>& vof = mf_vof.const_array(mfi);
        const amrex::Array4<const int>& iblank = mf_iblank.const_array(mfi);
        // Populate vof and density fluxes for each direction
        amrex::ParallelFor(
            xbx, ybx, zbx,
            [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                // vof flux
                amrex::Real f_vof =
                    alpha_flux(i, j, k, 0, i_th, vof, gpsi, iblank, dx);
                // check for overemptying depending on sign
                /*if (f_vof > 0.0) {
                    // vof is removed from cell i
                    f_vof = std::min(f_vof, vof(i, j, k) * dx[0] / dxf);
                } else {
                    // vof is removed from cell i-1
                    f_vof = std::max(f_vof, -vof(i - 1, j, k) * dx[0] / dxf);
                }*/
                // save modified flux
                fx(i, j, k, 0) = f_vof;
                // density flux
                // fx(i, j, k, 1) = rho1 * fx(i,j,k,0) + rho2 * (1.0 -
                // fx(i,j,k,0)); pressure gradient fluxes
            },
            [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                amrex::Real f_vof =
                    alpha_flux(i, j, k, 1, i_th, vof, gpsi, iblank, dx);
                // check for overemptying depending on sign
                /*if (f_vof > 0.0) {
                    // vof is removed from cell j
                    f_vof = std::min(f_vof, vof(i, j, k) * dx[1] / dxf);
                } else {
                    // vof is removed from cell j-1
                    f_vof = std::max(f_vof, -vof(i, j - 1, k) * dx[1] / dxf);
                }*/
                fy(i, j, k, 0) = f_vof;
            },
            [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                amrex::Real f_vof =
                    alpha_flux(i, j, k, 2, i_th, vof, gpsi, iblank, dx);
                /*if (i == 20 && j == 20 && k == 33)
                    std::cout << "f_vof " << f_vof;
                // check for overemptying depending on sign
                if (f_vof > 0.0) {
                    // vof is removed from cell k
                    f_vof = std::min(f_vof, vof(i, j, k) * dx[2] / dxf);
                } else {
                    // vof is removed from cell k-1
                    f_vof = std::max(f_vof, -vof(i, j, k - 1) * dx[2] / dxf);
                }*/
                if (i == 20 && j == 20 && k == 33)
                    std::cout << " " << f_vof << std::endl;
                fz(i, j, k, 0) = f_vof;
            });
    }
}

void process_fluxes(
    amrex::MultiFab& mf_fx,
    amrex::MultiFab& mf_fy,
    amrex::MultiFab& mf_fz,
    amrex::MultiFab& mf_vof,
    amrex::iMultiFab& mf_iblank)
{
    int ncomp = mf_fx.nComp();
    constexpr amrex::Real tiny = 1e-12;
    for (amrex::MFIter mfi(mf_iblank); mfi.isValid(); ++mfi) {
        const auto& vbx = mfi.validbox();
        const auto& xbx = amrex::surroundingNodes(vbx, 0);
        const auto& ybx = amrex::surroundingNodes(vbx, 1);
        const auto& zbx = amrex::surroundingNodes(vbx, 2);
        const amrex::Array4<amrex::Real>& fx = mf_fx.array(mfi);
        const amrex::Array4<amrex::Real>& fy = mf_fy.array(mfi);
        const amrex::Array4<amrex::Real>& fz = mf_fz.array(mfi);
        const amrex::Array4<const amrex::Real>& vof = mf_vof.const_array(mfi);
        const amrex::Array4<const int>& iblank = mf_iblank.const_array(mfi);
        // Only faces with iblank = -1 on both sides can have nonzero flux
        // Faces with single-phase on one side cannot have flux
        amrex::ParallelFor(
            xbx, ybx, zbx,
            [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                bool zero_all = (iblank(i - 1, j, k) + iblank(i, j, k) > -2);
                /*zero_all = zero_all ||
                           ((vof(i - 1, j, k) < tiny ||
                             vof(i - 1, j, k) > 1. - tiny) ||
                            (vof(i, j, k) < tiny || vof(i, j, k) > 1. -
                   tiny));*/
                for (int n = 0; n < ncomp; ++n) {
                    fx(i, j, k, n) *= zero_all ? 0. : 1.;
                }
            },
            [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                bool zero_all = (iblank(i, j - 1, k) + iblank(i, j, k) > -2);
                /*zero_all = zero_all ||
                           ((vof(i, j - 1, k) < tiny ||
                             vof(i, j - 1, k) > 1. - tiny) ||
                            (vof(i, j, k) < tiny || vof(i, j, k) > 1. -
                   tiny));*/
                for (int n = 0; n < ncomp; ++n) {
                    fy(i, j, k, n) *= zero_all ? 0. : 1.;
                }
            },
            [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                bool zero_all = (iblank(i, j, k - 1) + iblank(i, j, k) > -2);
                /*zero_all = zero_all ||
                           ((vof(i, j, k - 1) < tiny ||
                             vof(i, j, k - 1) > 1. - tiny) ||
                            (vof(i, j, k) < tiny || vof(i, j, k) > 1. -
                   tiny));*/
                for (int n = 0; n < ncomp; ++n) {
                    fz(i, j, k, n) *= zero_all ? 0. : 1.;
                }
            });
    }
}

void apply_fluxes(
    amrex::MultiFab& mf_fx,
    amrex::MultiFab& mf_fy,
    amrex::MultiFab& mf_fz,
    amrex::MultiFab& mf_vof,
    amrex::iMultiFab& mf_iblank,
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dx,
    const amrex::Real i_th)
{
    int ncomp = mf_fx.nComp();
    constexpr amrex::Real tiny = 1e-12;
    for (amrex::MFIter mfi(mf_vof); mfi.isValid(); ++mfi) {
        const auto& vbx = mfi.validbox();
        const amrex::Array4<amrex::Real>& fx = mf_fx.array(mfi);
        const amrex::Array4<amrex::Real>& fy = mf_fy.array(mfi);
        const amrex::Array4<amrex::Real>& fz = mf_fz.array(mfi);
        const amrex::Array4<amrex::Real>& vof = mf_vof.array(mfi);
        const amrex::Array4<const int>& iblank = mf_iblank.const_array(mfi);
        // Only faces with iblank = -1 on both sides can have nonzero flux
        amrex::ParallelFor(
            vbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                const amrex::Real oldvof = vof(i, j, k);
                vof(i, j, k) +=
                    i_th *
                    ((iblank(i, j, k) != -1)
                         ? 0.
                         : (fx(i + 1, j, k, 0) - fx(i, j, k, 0)) / dx[0] +
                               (fy(i, j + 1, k, 0) - fy(i, j, k, 0)) / dx[1] +
                               (fz(i, j, k + 1, 0) - fz(i, j, k, 0)) / dx[2]);
                // if (i == 20 && j == 20 && k == 33) {
                if (i == 20 && j == 20 && vof(i, j, k) > 1.0) {
                    std::cout << i << " " << j << " " << k << " "
                              << vof(i, j, k) << " " << oldvof << " "
                              << fz(i, j, k + 1, 0) << " " << fz(i, j, k, 0)
                              << std::endl;
                }
                /*if (vof(i, j, k) > 0.0 && vof(i, j, k) < 1e-8)
                    vof(i, j, k) = 0.0;*/
                // vof(i, j, k) = std::min(1.0, std::max(0.0, vof(i, j, k)));

                // This outer conditional is just for the sake of debugging
                if (vof(i, j, k) < 1.0 && vof(i, j, k) > 0.0) {
                    vof(i, j, k) =
                        vof(i, j, k) < tiny
                            ? 0.0
                            : (vof(i, j, k) > 1. - tiny ? 1. : vof(i, j, k));
                }
            });
    }
}

// for debugging
void equate_vof(amrex::MultiFab& mf_gpsi, amrex::MultiFab& mf_vof)
{
    for (amrex::MFIter mfi(mf_vof); mfi.isValid(); ++mfi) {
        const auto& vbx = mfi.validbox();
        const amrex::Array4<const amrex::Real>& gpsi = mf_gpsi.const_array(mfi);
        const amrex::Array4<amrex::Real>& vof = mf_vof.array(mfi);
        // Only faces with iblank = -1 on both sides can have nonzero flux
        amrex::ParallelFor(
            vbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                vof(i, j, k) = gpsi(i, j, k, 3);
            });
    }
}

void SharpenNaluData(CFDSim& sim, const int nit)
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

    // Get blanking for cells
    auto& iblank_cell = repo.get_int_field("iblank_cell");

    // Get fields that will be modified
    auto& vof = repo.get_field("vof");
    auto& rho = repo.get_field("density");
    auto& gp = repo.get_field("gp");

    // Create scratch fields for fluxes
    // (components are vof, density, gpx, gpy, gpz)
    auto flux_x = repo.create_scratch_field(5, 0, amr_wind::FieldLoc::XFACE);
    auto flux_y = repo.create_scratch_field(5, 0, amr_wind::FieldLoc::YFACE);
    auto flux_z = repo.create_scratch_field(5, 0, amr_wind::FieldLoc::ZFACE);
    // Create scratch field for approximate signed distance function and grad
    // (components 0-2 are gradient, 3 is asdf)
    auto grad_asdf = repo.create_scratch_field(4, vof.num_grow()[0]);

    // Process vof for tiny margins from single-phase
    for (int lev = 0; lev < nlevels; ++lev) {
        process_vof(vof(lev));
    }

    // Pseudo-time loop
    for (int n = 0; n < nit; ++n) {
        std::cout << "new step " << n << std::endl;
        for (int lev = 0; lev < nlevels; ++lev) {
            // Thickness of interface should be 1/2 the size of the cell?
            auto dx = (geom[lev]).CellSizeArray();
            const amrex::Real i_th = 0.5 * std::cbrt(dx[0] * dx[1] * dx[2]);

            // Populate approximate signed distance function
            populate_grad_psi(
                (*grad_asdf)(lev), vof(lev), iblank_cell(lev), dx, i_th);

            // Sharpening fluxes for vof, density, pressure gradient
            populate_sharpen_fluxes(
                (*flux_x)(lev), (*flux_y)(lev), (*flux_z)(lev),
                (*grad_asdf)(lev), vof(lev), iblank_cell(lev), dx, i_th, rhol,
                rhog);

            // Zero fluxes outside nalu portion of domain
            process_fluxes(
                (*flux_x)(lev), (*flux_y)(lev), (*flux_z)(lev), vof(lev),
                iblank_cell(lev));

            // Apply fluxes
            apply_fluxes(
                (*flux_x)(lev), (*flux_y)(lev), (*flux_z)(lev), vof(lev),
                iblank_cell(lev), dx, i_th);

            // Replace vof with another field to look at it
            // equate_vof((*grad_asdf)(lev), vof(lev));
        }

        // Average down fluxes
        // Get pseudo CFL (multiplier for sharpening term)
        // Apply fluxes

        // Average down vof for consistency

        // Fillpatch for ghost cells
        vof.fillpatch(sim.time().current_time());
    }
}
} // namespace overset
} // namespace amr_wind

#endif