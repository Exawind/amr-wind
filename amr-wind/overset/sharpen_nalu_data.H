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

void process_vof(amrex::MultiFab& mf_vof, const amrex::Real vof_tol)
{
    for (amrex::MFIter mfi(mf_vof); mfi.isValid(); ++mfi) {
        const auto& gbx = mfi.growntilebox();
        const amrex::Array4<amrex::Real>& vof = mf_vof.array(mfi);
        // Populate approximate signed distance function
        amrex::ParallelFor(
            gbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                // This conditional is just for the sake of debugging
                if (vof(i, j, k) < 1.0 && vof(i, j, k) > 0.0) {
                    vof(i, j, k) =
                        vof(i, j, k) < vof_tol
                            ? 0.0
                            : (vof(i, j, k) > 1. - vof_tol ? 1. : vof(i, j, k));
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

amrex::Real alpha_discrete_flux(
    int i,
    int j,
    int k,
    int dir,
    const amrex::Real margin,
    amrex::Array4<amrex::Real const> const& vof,
    amrex::Array4<amrex::Real const> const& tg_vof,
    amrex::Array4<amrex::Real const> const& normal)
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
}

void populate_sharpen_discrete_fluxes(
    amrex::MultiFab& mf_fx,
    amrex::MultiFab& mf_fy,
    amrex::MultiFab& mf_fz,
    amrex::MultiFab& mf_vof,
    amrex::MultiFab& mf_target_vof,
    amrex::MultiFab& mf_norm,
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
        const amrex::Array4<const amrex::Real>& vel =
            mf_velocity.const_array(mfi);
        // Populate vof and density fluxes for each direction
        amrex::ParallelFor(
            xbx, ybx, zbxp1,
            [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                // vof flux
                amrex::Real flux =
                    alpha_discrete_flux(i, j, k, 0, margin, vof, tg_vof, norm);
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
                amrex::Real flux =
                    alpha_discrete_flux(i, j, k, 1, margin, vof, tg_vof, norm);
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
                amrex::Real flux =
                    alpha_discrete_flux(i, j, k, 2, margin, vof, tg_vof, norm);
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
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dx)
{
    int ncomp = mf_fx.nComp();
    constexpr amrex::Real tiny = 1e-12;
    // Pseudo-time factor
    const amrex::Real ptfac = 1.0;
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
                vof(i, j, k) += ptfac * ((fx(i + 1, j, k, 0) - fx(i, j, k, 0)) +
                                         (fy(i, j + 1, k, 0) - fy(i, j, k, 0)) +
                                         (fz(i, j, k + 1, 0) - fz(i, j, k, 0)));
                dens(i, j, k) +=
                    ptfac * ((fx(i + 1, j, k, 1) - fx(i, j, k, 1)) +
                             (fy(i, j + 1, k, 1) - fy(i, j, k, 1)) +
                             (fz(i, j, k + 1, 1) - fz(i, j, k, 1)));
                vel(i, j, k, 0) =
                    1.0 / dens(i, j, k) *
                    (olddens * vel(i, j, k, 0) +
                     ptfac * ((fx(i + 1, j, k, 2) - fx(i, j, k, 2)) +
                              (fy(i, j + 1, k, 2) - fy(i, j, k, 2)) +
                              (fz(i, j, k + 1, 2) - fz(i, j, k, 2))));
                vel(i, j, k, 1) =
                    1.0 / dens(i, j, k) *
                    (olddens * vel(i, j, k, 1) +
                     ptfac * ((fx(i + 1, j, k, 3) - fx(i, j, k, 3)) +
                              (fy(i, j + 1, k, 3) - fy(i, j, k, 3)) +
                              (fz(i, j, k + 1, 3) - fz(i, j, k, 3))));
                vel(i, j, k, 2) =
                    1.0 / dens(i, j, k) *
                    (olddens * vel(i, j, k, 2) +
                     ptfac * ((fx(i + 1, j, k, 4) - fx(i, j, k, 4)) +
                              (fy(i, j + 1, k, 4) - fy(i, j, k, 4)) +
                              (fz(i, j, k + 1, 4) - fz(i, j, k, 4))));

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
                p(i, j, k) += ptfac * sp(i, j, k) * dx[2];
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

void replace_gradp_hs(
    amrex::MultiFab& mf_gp,
    amrex::MultiFab& mf_density,
    amrex::iMultiFab& mf_iblank,
    const amrex::Real grav_z)
{
    for (amrex::MFIter mfi(mf_gp); mfi.isValid(); ++mfi) {
        const auto& gbx = mfi.growntilebox();
        const amrex::Array4<amrex::Real>& gp = mf_gp.array(mfi);
        const amrex::Array4<const amrex::Real>& rho =
            mf_density.const_array(mfi);
        const amrex::Array4<const int>& iblank = mf_iblank.const_array(mfi);
        amrex::ParallelFor(
            gbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                // Replace amr-wind vof values with originals
                if (iblank(i, j, k) == -1) {
                    gp(i, j, k, 0) = 0.;
                    gp(i, j, k, 1) = 0.;
                    gp(i, j, k, 2) = rho(i, j, k) * grav_z;
                }
            });
    }
}

void SharpenNaluDataDiscrete(
    CFDSim& sim,
    const int nit,
    const amrex::Real tol,
    const int ncconv,
    const amrex::Real rlscale,
    const amrex::Real margin,
    const amrex::Real proc_tgvof_tol,
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

    // Process target vof for tiny margins from single-phase
    for (int lev = 0; lev < nlevels; ++lev) {
        // A tolerance of 0 should do nothing
        process_vof((*target_vof)(lev), proc_tgvof_tol);
    }

    // Replace vof with original values in amr domain
    for (int lev = 0; lev < nlevels; ++lev) {
        harmonize_vof((*target_vof)(lev), vof(lev), iblank_cell(lev));
    }

    // Purely for debugging, should be removed later
    for (int lev = 0; lev < nlevels; ++lev) {
        equate_vof((*target_vof)(lev), levelset(lev));
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
            // Height of cell is needed for pressure source
            auto dx = (geom[lev]).CellSizeArray();

            // Populate normal vector
            populate_normal_vector(
                (*normal_vec)(lev), vof(lev), iblank_cell(lev));

            // Sharpening fluxes for vof, density, and momentum
            populate_sharpen_discrete_fluxes(
                (*flux_x)(lev), (*flux_y)(lev), (*flux_z)(lev), vof(lev),
                (*target_vof)(lev), (*normal_vec)(lev), velocity(lev), margin,
                rhol, rhog);

            // Process fluxes and get pressure source term
            process_fluxes_calc_src(
                (*flux_x)(lev), (*flux_y)(lev), (*flux_z)(lev), (*p_src)(lev),
                iblank_cell(lev), grav_z);

            // Apply fluxes
            apply_fluxes(
                (*flux_x)(lev), (*flux_y)(lev), (*flux_z)(lev), (*p_src)(lev),
                vof(lev), rho(lev), velocity(lev), p(lev), dx);

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

void ReplaceMaskedGradP(CFDSim& sim)
{
    auto& repo = sim.repo();
    auto nlevels = repo.num_active_levels();
    auto geom = sim.mesh().Geom();

    // Get gravity
    auto& mphase = sim.physics_manager().get<MultiPhase>();
    const amrex::Real grav_z = mphase.gravity()[2];

    // Get blanking for cells
    auto& iblank_cell = repo.get_int_field("iblank_cell");

    // Get fields that will be modified or used
    auto& rho = repo.get_field("density");
    auto& gp = repo.get_field("gp");

    // Replace initial gp with best guess (hydrostatic)
    for (int lev = 0; lev < nlevels; ++lev) {
        // Populate approximate signed distance function
        replace_gradp_hs(gp(lev), rho(lev), iblank_cell(lev), grav_z);
    }
}
} // namespace overset
} // namespace amr_wind

#endif