#ifndef SHARPEN_NALU_DATA_H_
#define SHARPEN_NALU_DATA_H_

#include "amr-wind/CFDSim.H"
#include "amr-wind/physics/multiphase/MultiPhase.H"
#include "amr-wind/equation_systems/vof/volume_fractions.H"
#include "amr-wind/utilities/IOManager.H"

namespace amr_wind::overset {

// Approximate signed distance function
amrex::Real AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
asdf(const amrex::Real a_vof, const amrex::Real i_th)
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
        amrex::ParallelFor(
            gbxm1, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                // Neumann condition across nalu bdy
                int ibdy = 0;
                int jbdy = 0;
                int kbdy = 0;
                if (iblank(i, j, k) != iblank(i - 1, j, k)) {
                    ibdy = -1;
                }
                if (iblank(i, j, k) != iblank(i + 1, j, k)) {
                    ibdy += 1;
                }
                if (iblank(i, j, k) != iblank(i, j - 1, k)) {
                    jbdy = -1;
                }
                if (iblank(i, j, k) != iblank(i, j + 1, k)) {
                    jbdy += 1;
                }
                // Assuming gravity is in z
                /*if (iblank(i, j, k) != iblank(i, j, k - 1)) {
                    kbdy = -1;
                }
                if (iblank(i, j, k) != iblank(i, j, k + 1)) {
                    kbdy += 1;
                }*/
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

amrex::Real AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE alpha_discrete_flux(
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

amrex::Real AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE alpha_discrete_flux_bndry(
    int i,
    int j,
    int k,
    int dir,
    amrex::Array4<amrex::Real const> const& vof,
    amrex::Array4<amrex::Real const> const& tg_vof,
    amrex::Array4<amrex::Real const> const& normal,
    amrex::Array4<int const> const& iblank)
{
    // Set up neighbor indices
    int ii = i;
    int jj = j;
    int kk = k;
    ii += (dir == 0) ? -1 : 0;
    jj += (dir == 1) ? -1 : 0;
    kk += (dir == 2) ? -1 : 0;

    // Normal vector in each cell (already normalized)
    const amrex::Real norm_ = normal(i, j, k, dir);
    const amrex::Real norm_nb = normal(ii, jj, kk, dir);

    // Determine which delta_phi (and multiply by normal)
    // The sign depends on side of flux face (like upwinding)
    const amrex::Real dphi_ = (tg_vof(i, j, k) - vof(i, j, k)) * (-norm_);
    const amrex::Real dphi_nb =
        (tg_vof(ii, jj, kk) - vof(ii, jj, kk)) * norm_nb;
    // Average value for initialization
    amrex::Real dphi_eval = 0.5 * (dphi_ + dphi_nb);
    // Always prefer the dphi where iblank = -1
    // Reduce the magnitude of these fluxes to limit consequence on domain
    constexpr amrex::Real reduce_factor = 0.5;
    if (iblank(i, j, k) == -1 && iblank(ii, jj, kk) == 1) {
        dphi_eval = reduce_factor * dphi_;
    } else if (iblank(i, j, k) == 1 && iblank(ii, jj, kk) == -1) {
        dphi_eval = reduce_factor * dphi_nb;
    }
    return dphi_eval;
}

void AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE velocity_discrete_face(
    int i,
    int j,
    int k,
    int dir,
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
    // Average value when gphi = 0
    uface = 0.5 * (u_ + u_nb);
    vface = 0.5 * (v_ + v_nb);
    wface = 0.5 * (w_ + w_nb);
    // Use simple upwinding
    if (gphi > 0.0) {
        uface = u_nb;
        vface = v_nb;
        wface = w_nb;
    }
    if (gphi < 0.0) {
        uface = u_;
        vface = v_;
        wface = w_;
    }
}

void AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE velocity_discrete_face_bndry(
    int i,
    int j,
    int k,
    int dir,
    amrex::Array4<int const> const& iblank,
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

    // Get velocities on both sides
    const amrex::Real u_ = velocity(i, j, k, 0);
    const amrex::Real v_ = velocity(i, j, k, 1);
    const amrex::Real w_ = velocity(i, j, k, 2);
    const amrex::Real u_nb = velocity(ii, jj, kk, 0);
    const amrex::Real v_nb = velocity(ii, jj, kk, 1);
    const amrex::Real w_nb = velocity(ii, jj, kk, 2);
    // Average value for initialization
    uface = 0.5 * (u_ + u_nb);
    vface = 0.5 * (v_ + v_nb);
    wface = 0.5 * (w_ + w_nb);
    // Defer to iblank == -1 values
    if (iblank(i, j, k) == -1 && iblank(ii, jj, kk) == 1) {
        uface = u_;
        vface = v_;
        wface = w_;
    } else if (iblank(i, j, k) == 1 && iblank(ii, jj, kk) == -1) {
        uface = u_nb;
        vface = v_nb;
        wface = w_nb;
    }
}

void AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE gp_rho_discrete_face(
    int i,
    int j,
    int k,
    int dir,
    amrex::Array4<amrex::Real const> const& vof,
    amrex::Array4<amrex::Real const> const& gp,
    amrex::Array4<amrex::Real const> const& rho,
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
    const amrex::Real u_ = gp(i, j, k, 0) / rho(i, j, k);
    const amrex::Real v_ = gp(i, j, k, 1) / rho(i, j, k);
    const amrex::Real w_ = gp(i, j, k, 2) / rho(i, j, k);
    const amrex::Real u_nb = gp(ii, jj, kk, 0) / rho(ii, jj, kk);
    const amrex::Real v_nb = gp(ii, jj, kk, 1) / rho(ii, jj, kk);
    const amrex::Real w_nb = gp(ii, jj, kk, 2) / rho(ii, jj, kk);
    // Average value when gphi = 0
    uface = 0.5 * (u_ + u_nb);
    vface = 0.5 * (v_ + v_nb);
    wface = 0.5 * (w_ + w_nb);
    // Use simple upwinding
    if (gphi > 0.0) {
        uface = u_nb;
        vface = v_nb;
        wface = w_nb;
    }
    if (gphi < 0.0) {
        uface = u_;
        vface = v_;
        wface = w_;
    }
}

void AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE gp_rho_discrete_face_bndry(
    int i,
    int j,
    int k,
    int dir,
    amrex::Array4<int const> const& iblank,
    amrex::Array4<amrex::Real const> const& gp,
    amrex::Array4<amrex::Real const> const& rho,
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

    // Get velocities on both sides
    const amrex::Real u_ = gp(i, j, k, 0) / rho(i, j, k);
    const amrex::Real v_ = gp(i, j, k, 1) / rho(i, j, k);
    const amrex::Real w_ = gp(i, j, k, 2) / rho(i, j, k);
    const amrex::Real u_nb = gp(ii, jj, kk, 0) / rho(ii, jj, kk);
    const amrex::Real v_nb = gp(ii, jj, kk, 1) / rho(ii, jj, kk);
    const amrex::Real w_nb = gp(ii, jj, kk, 2) / rho(ii, jj, kk);
    // Average value for initialization
    uface = 0.5 * (u_ + u_nb);
    vface = 0.5 * (v_ + v_nb);
    wface = 0.5 * (w_ + w_nb);
    // Defer to iblank == -1 values
    if (iblank(i, j, k) == -1 && iblank(ii, jj, kk) == 1) {
        uface = u_;
        vface = v_;
        wface = w_;
    } else if (iblank(i, j, k) == 1 && iblank(ii, jj, kk) == -1) {
        uface = u_nb;
        vface = v_nb;
        wface = w_nb;
    }
}

void populate_sharpen_discrete_fluxes(
    amrex::iMultiFab& mf_iblank,
    amrex::MultiFab& mf_fx,
    amrex::MultiFab& mf_fy,
    amrex::MultiFab& mf_fz,
    amrex::MultiFab& mf_vof,
    amrex::MultiFab& mf_target_vof,
    amrex::MultiFab& mf_norm,
    amrex::MultiFab& mf_velocity,
    amrex::MultiFab& mf_gp,
    amrex::MultiFab& mf_density,
    const amrex::Real Gamma,
    const amrex::Real margin,
    const amrex::Real rho1,
    const amrex::Real rho2)
{
    for (amrex::MFIter mfi(mf_vof); mfi.isValid(); ++mfi) {
        const auto& vbx = mfi.validbox();
        const auto& xbxp1 = grow(amrex::surroundingNodes(vbx, 0), 1);
        const auto& ybxp1 = grow(amrex::surroundingNodes(vbx, 1), 1);
        const auto& zbxp1 = grow(amrex::surroundingNodes(vbx, 2), 1);
        const amrex::Array4<int>& iblank = mf_iblank.array(mfi);
        const amrex::Array4<amrex::Real>& fx = mf_fx.array(mfi);
        const amrex::Array4<amrex::Real>& fy = mf_fy.array(mfi);
        const amrex::Array4<amrex::Real>& fz = mf_fz.array(mfi);
        const amrex::Array4<const amrex::Real>& vof = mf_vof.const_array(mfi);
        const amrex::Array4<const amrex::Real>& tg_vof =
            mf_target_vof.const_array(mfi);
        const amrex::Array4<const amrex::Real>& norm = mf_norm.const_array(mfi);
        const amrex::Array4<const amrex::Real>& vel =
            mf_velocity.const_array(mfi);
        /*const amrex::Array4<const amrex::Real>& gp_rho =
            mf_gp_rho.const_array(mfi);*/
        const amrex::Array4<const amrex::Real>& gp = mf_gp.const_array(mfi);
        const amrex::Array4<const amrex::Real>& rho =
            mf_density.const_array(mfi);
        // Populate vof and density fluxes for each direction
        amrex::ParallelFor(
            xbxp1, ybxp1, zbxp1,
            [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                // vof flux
                amrex::Real flux =
                    Gamma *
                    alpha_discrete_flux(i, j, k, 0, margin, vof, tg_vof, norm);
                fx(i, j, k, 0) = flux;
                // density flux
                flux *= (rho1 - rho2);
                fx(i, j, k, 1) = flux;
                // momentum fluxes (dens flux * face vel)
                amrex::Real uf, vf, wf;
                velocity_discrete_face(i, j, k, 0, vof, vel, uf, vf, wf);
                fx(i, j, k, 2) = flux * uf;
                fx(i, j, k, 3) = flux * vf;
                fx(i, j, k, 4) = flux * wf;
                gp_rho_discrete_face(i, j, k, 0, vof, gp, rho, uf, vf, wf);
                fx(i, j, k, 5) = flux * uf;
                fx(i, j, k, 6) = flux * vf;
                fx(i, j, k, 7) = flux * wf;
                // Turn "on" all flux faces, later modified in process_fluxes
                fx(i, j, k, 8) = 1.0;
            },
            [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                amrex::Real flux =
                    Gamma *
                    alpha_discrete_flux(i, j, k, 1, margin, vof, tg_vof, norm);
                fy(i, j, k, 0) = flux;
                flux *= (rho1 - rho2);
                fy(i, j, k, 1) = flux;
                amrex::Real uf, vf, wf;
                velocity_discrete_face(i, j, k, 1, vof, vel, uf, vf, wf);
                fy(i, j, k, 2) = flux * uf;
                fy(i, j, k, 3) = flux * vf;
                fy(i, j, k, 4) = flux * wf;
                gp_rho_discrete_face(i, j, k, 1, vof, gp, rho, uf, vf, wf);
                fy(i, j, k, 5) = flux * uf;
                fy(i, j, k, 6) = flux * vf;
                fy(i, j, k, 7) = flux * wf;
                fy(i, j, k, 8) = 1.0;
            },
            [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                amrex::Real flux =
                    Gamma *
                    alpha_discrete_flux(i, j, k, 2, margin, vof, tg_vof, norm);
                if (iblank(i, j, k) + iblank(i, j, k - 1) > -2) {
                    // Calculate with modified flux value across overset bndry
                    flux = Gamma * alpha_discrete_flux_bndry(
                                       i, j, k, 2, vof, tg_vof, norm, iblank);
                }
                fz(i, j, k, 0) = flux;
                flux *= (rho1 - rho2);
                fz(i, j, k, 1) = flux;
                amrex::Real uf, vf, wf;
                velocity_discrete_face(i, j, k, 2, vof, vel, uf, vf, wf);
                if (iblank(i, j, k) + iblank(i, j, k - 1) > -2) {
                    // Calculate with modified flux value across overset bndry
                    velocity_discrete_face_bndry(
                        i, j, k, 2, iblank, vel, uf, vf, wf);
                }
                fz(i, j, k, 2) = flux * uf;
                fz(i, j, k, 3) = flux * vf;
                fz(i, j, k, 4) = flux * wf;
                gp_rho_discrete_face(i, j, k, 2, vof, gp, rho, uf, vf, wf);
                if (iblank(i, j, k) + iblank(i, j, k - 1) > -2) {
                    // Calculate with modified flux value across overset bndry
                    gp_rho_discrete_face_bndry(
                        i, j, k, 2, iblank, gp, rho, uf, vf, wf);
                }
                fz(i, j, k, 5) = flux * uf;
                fz(i, j, k, 6) = flux * vf;
                fz(i, j, k, 7) = flux * wf;
                fz(i, j, k, 8) = 1.0;
            });
    }
}

void process_fluxes_calc_src(
    amrex::MultiFab& mf_fx,
    amrex::MultiFab& mf_fy,
    amrex::MultiFab& mf_fz,
    amrex::MultiFab& mf_vof,
    amrex::MultiFab& mf_psource,
    amrex::iMultiFab& mf_iblank,
    const bool sh_p_flag)
{
    int ncompx = mf_fx.nComp();
    int ncompy = mf_fy.nComp();
    int ncompz = mf_fz.nComp();
    if (!sh_p_flag) {
        mf_psource.setVal(0.0);
    }
    for (amrex::MFIter mfi(mf_iblank); mfi.isValid(); ++mfi) {
        const auto& vbx = mfi.validbox();
        const auto& xbxp1 = grow(amrex::surroundingNodes(vbx, 0), 1);
        const auto& ybxp1 = grow(amrex::surroundingNodes(vbx, 1), 1);
        const auto& zbxp1 = grow(amrex::surroundingNodes(vbx, 2), 1);
        const amrex::Array4<amrex::Real>& fx = mf_fx.array(mfi);
        const amrex::Array4<amrex::Real>& fy = mf_fy.array(mfi);
        const amrex::Array4<amrex::Real>& fz = mf_fz.array(mfi);
        const amrex::Array4<amrex::Real>& vof = mf_vof.array(mfi);
        const amrex::Array4<amrex::Real>& sp = mf_psource.array(mfi);
        const amrex::Array4<const int>& iblank = mf_iblank.const_array(mfi);
        // Only faces with iblank = -1 on both sides can have nonzero flux
        amrex::ParallelFor(
            xbxp1, ybxp1, zbxp1,
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
                // Fluxing is allowed across overset boundaries in z
                bool zero_all =
                    std::min(iblank(i, j, k - 1), iblank(i, j, k)) > -1;
                for (int n = 0; n < ncompz; ++n) {
                    fz(i, j, k, n) *= zero_all ? 0. : 1.;
                }
            });
        // With knowledge of fluxes, compute pressure source term
        amrex::Box const& nbx = mfi.nodaltilebox();
        if (sh_p_flag) {
            constexpr amrex::Real tiny = 1e-100;
            amrex::ParallelFor(
                nbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    // Get components of outer product between f and gradp/rho
                    // from interpolating between values at neighboring faces
                    const amrex::Real fxgpx_flux =
                        (fx(i, j, k, 5) + fx(i, j, k - 1, 5) +
                         fx(i, j - 1, k, 5) + fx(i, j - 1, k - 1, 5)) /
                        (fx(i, j, k, 8) + fx(i, j, k - 1, 8) +
                         fx(i, j - 1, k, 8) + fx(i, j - 1, k - 1, 8) + tiny);
                    const amrex::Real fxgpy_flux =
                        (fx(i, j, k, 6) + fx(i, j, k - 1, 6) +
                         fx(i, j - 1, k, 6) + fx(i, j - 1, k - 1, 6)) /
                        (fx(i, j, k, 8) + fx(i, j, k - 1, 8) +
                         fx(i, j - 1, k, 8) + fx(i, j - 1, k - 1, 8) + tiny);
                    const amrex::Real fxgpz_flux =
                        (fx(i, j, k, 7) + fx(i, j, k - 1, 7) +
                         fx(i, j - 1, k, 7) + fx(i, j - 1, k - 1, 7)) /
                        (fx(i, j, k, 8) + fx(i, j, k - 1, 8) +
                         fx(i, j - 1, k, 8) + fx(i, j - 1, k - 1, 8) + tiny);

                    const amrex::Real fygpx_flux =
                        (fy(i, j, k, 5) + fy(i - 1, j, k, 5) +
                         fy(i, j, k - 1, 5) + fy(i - 1, j, k - 1, 5)) /
                        (fy(i, j, k, 8) + fy(i - 1, j, k, 8) +
                         fy(i, j, k - 1, 8) + fy(i - 1, j, k - 1, 8) + tiny);
                    const amrex::Real fygpy_flux =
                        (fy(i, j, k, 6) + fy(i - 1, j, k, 6) +
                         fy(i, j, k - 1, 6) + fy(i - 1, j, k - 1, 6)) /
                        (fy(i, j, k, 8) + fy(i - 1, j, k, 8) +
                         fy(i, j, k - 1, 8) + fy(i - 1, j, k - 1, 8) + tiny);
                    const amrex::Real fygpz_flux =
                        (fy(i, j, k, 7) + fy(i - 1, j, k, 7) +
                         fy(i, j, k - 1, 7) + fy(i - 1, j, k - 1, 7)) /
                        (fy(i, j, k, 8) + fy(i - 1, j, k, 8) +
                         fy(i, j, k - 1, 8) + fy(i - 1, j, k - 1, 8) + tiny);

                    const amrex::Real fzgpx_flux =
                        (fz(i, j, k, 5) + fz(i - 1, j, k, 5) +
                         fz(i, j - 1, k, 5) + fz(i - 1, j - 1, k, 5)) /
                        (fz(i, j, k, 8) + fz(i - 1, j, k, 8) +
                         fz(i, j - 1, k, 8) + fz(i - 1, j - 1, k, 8) + tiny);
                    const amrex::Real fzgpy_flux =
                        (fz(i, j, k, 6) + fz(i - 1, j, k, 6) +
                         fz(i, j - 1, k, 6) + fz(i - 1, j - 1, k, 6)) /
                        (fz(i, j, k, 8) + fz(i - 1, j, k, 8) +
                         fz(i, j - 1, k, 8) + fz(i - 1, j - 1, k, 8) + tiny);
                    const amrex::Real fzgpz_flux =
                        (fz(i, j, k, 7) + fz(i - 1, j, k, 7) +
                         fz(i, j - 1, k, 7) + fz(i - 1, j - 1, k, 7)) /
                        (fz(i, j, k, 8) + fz(i - 1, j, k, 8) +
                         fz(i, j - 1, k, 8) + fz(i - 1, j - 1, k, 8) + tiny);

                    // Calculate interface normal at pressure node
                    amrex::Real n0 =
                        ((vof(i, j, k) - vof(i - 1, j, k)) * fx(i, j, k, 8) +
                         (vof(i, j - 1, k) - vof(i - 1, j - 1, k)) *
                             fx(i, j - 1, k, 8) +
                         (vof(i, j, k - 1) - vof(i - 1, j, k - 1)) *
                             fx(i, j, k - 1, 8) +
                         (vof(i, j - 1, k - 1) - vof(i - 1, j - 1, k - 1)) *
                             fx(i, j - 1, k - 1, 8)) /
                        (fx(i, j, k, 8) + fx(i, j, k - 1, 8) +
                         fx(i, j - 1, k, 8) + fx(i, j - 1, k - 1, 8) + tiny);
                    amrex::Real n1 =
                        ((vof(i, j, k) - vof(i, j - 1, k)) * fy(i, j, k, 8) +
                         (vof(i - 1, j, k) - vof(i - 1, j - 1, k)) *
                             fy(i - 1, j, k, 8) +
                         (vof(i, j, k - 1) - vof(i, j - 1, k - 1)) *
                             fy(i, j, k - 1, 8) +
                         (vof(i - 1, j, k - 1) - vof(i - 1, j - 1, k - 1)) *
                             fy(i - 1, j, k - 1, 8)) /
                        (fy(i, j, k, 8) + fy(i, j, k - 1, 8) +
                         fy(i - 1, j, k, 8) + fy(i - 1, j, k - 1, 8) + tiny);
                    amrex::Real n2 =
                        ((vof(i, j, k) - vof(i, j, k - 1)) * fz(i, j, k, 8) +
                         (vof(i - 1, j, k) - vof(i - 1, j, k - 1)) *
                             fz(i - 1, j, k, 8) +
                         (vof(i, j - 1, k) - vof(i, j - 1, k - 1)) *
                             fz(i, j - 1, k, 8) +
                         (vof(i - 1, j - 1, k) - vof(i - 1, j - 1, k - 1)) *
                             fz(i - 1, j - 1, k, 8)) /
                        (fz(i, j, k, 8) + fz(i, j - 1, k, 8) +
                         fz(i - 1, j, k, 8) + fz(i - 1, j - 1, k, 8) + tiny);
                    // Normalize vector
                    amrex::Real nmag =
                        std::sqrt(n0 * n0 + n1 * n1 + n2 * n2) + tiny;
                    n0 /= nmag;
                    n1 /= nmag;
                    n2 /= nmag;

                    // Perform double contraction
                    sp(i, j, k) = fxgpx_flux * n0 * n0 + fxgpy_flux * n0 * n1 +
                                  fxgpz_flux * n0 * n2 + fygpx_flux * n1 * n0 +
                                  fygpy_flux * n1 * n1 + fygpz_flux * n1 * n2 +
                                  fzgpx_flux * n2 * n0 + fzgpy_flux * n2 * n1 +
                                  fzgpz_flux * n2 * n2;
                });
        }
    }
}

void apply_fluxes(
    amrex::iMultiFab& mf_iblank,
    amrex::iMultiFab& mf_iblnod,
    amrex::MultiFab& mf_fx,
    amrex::MultiFab& mf_fy,
    amrex::MultiFab& mf_fz,
    amrex::MultiFab& mf_psource,
    amrex::MultiFab& mf_vof,
    amrex::MultiFab& mf_dens,
    amrex::MultiFab& mf_vel,
    amrex::MultiFab& mf_gp,
    amrex::MultiFab& mf_pressure,
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dx,
    amrex::Real ptfac,
    const bool mod_gp)
{
    constexpr amrex::Real tiny = 1e-12;
    // Pseudo-time factor
    // const amrex::Real ptfac = 1.0;
    for (amrex::MFIter mfi(mf_vof); mfi.isValid(); ++mfi) {
        const auto& vbx = mfi.validbox();
        amrex::Box const& nbx = mfi.nodaltilebox();
        const amrex::Array4<amrex::Real>& fx = mf_fx.array(mfi);
        const amrex::Array4<amrex::Real>& fy = mf_fy.array(mfi);
        const amrex::Array4<amrex::Real>& fz = mf_fz.array(mfi);
        const amrex::Array4<amrex::Real>& vof = mf_vof.array(mfi);
        const amrex::Array4<amrex::Real>& dens = mf_dens.array(mfi);
        const amrex::Array4<amrex::Real>& vel = mf_vel.array(mfi);
        const amrex::Array4<amrex::Real>& gp = mf_gp.array(mfi);
        const amrex::Array4<amrex::Real>& p = mf_pressure.array(mfi);
        const amrex::Array4<amrex::Real>& sp = mf_psource.array(mfi);
        const amrex::Array4<int>& iblank = mf_iblank.array(mfi);
        const amrex::Array4<int>& iblnod = mf_iblnod.array(mfi);
        amrex::ParallelFor(
            vbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                const amrex::Real olddens = dens(i, j, k);
                const amrex::Real chng = (iblank(i, j, k) == -1) ? 1.0 : 0.0;
                vof(i, j, k) += chng * ptfac *
                                ((fx(i + 1, j, k, 0) - fx(i, j, k, 0)) / dx[0] +
                                 (fy(i, j + 1, k, 0) - fy(i, j, k, 0)) / dx[1] +
                                 (fz(i, j, k + 1, 0) - fz(i, j, k, 0)) / dx[2]);
                dens(i, j, k) +=
                    chng * ptfac *
                    ((fx(i + 1, j, k, 1) - fx(i, j, k, 1)) / dx[0] +
                     (fy(i, j + 1, k, 1) - fy(i, j, k, 1)) / dx[1] +
                     (fz(i, j, k + 1, 1) - fz(i, j, k, 1)) / dx[2]);
                vel(i, j, k, 0) =
                    (1.0 - chng) * vel(i, j, k, 0) +
                    chng * 1.0 / dens(i, j, k) *
                        (olddens * vel(i, j, k, 0) +
                         ptfac *
                             ((fx(i + 1, j, k, 2) - fx(i, j, k, 2)) / dx[0] +
                              (fy(i, j + 1, k, 2) - fy(i, j, k, 2)) / dx[1] +
                              (fz(i, j, k + 1, 2) - fz(i, j, k, 2)) / dx[2]));
                vel(i, j, k, 1) =
                    (1.0 - chng) * vel(i, j, k, 0) +
                    chng * 1.0 / dens(i, j, k) *
                        (olddens * vel(i, j, k, 1) +
                         ptfac *
                             ((fx(i + 1, j, k, 3) - fx(i, j, k, 3)) / dx[0] +
                              (fy(i, j + 1, k, 3) - fy(i, j, k, 3)) / dx[1] +
                              (fz(i, j, k + 1, 3) - fz(i, j, k, 3)) / dx[2]));
                vel(i, j, k, 2) =
                    (1.0 - chng) * vel(i, j, k, 0) +
                    chng * 1.0 / dens(i, j, k) *
                        (olddens * vel(i, j, k, 2) +
                         ptfac *
                             ((fx(i + 1, j, k, 4) - fx(i, j, k, 4)) / dx[0] +
                              (fy(i, j + 1, k, 4) - fy(i, j, k, 4)) / dx[1] +
                              (fz(i, j, k + 1, 4) - fz(i, j, k, 4)) / dx[2]));

                if (mod_gp) {

                    gp(i, j, k, 0) +=
                        chng * ptfac *
                        ((fx(i + 1, j, k, 5) - fx(i, j, k, 5)) / dx[0] +
                         (fy(i, j + 1, k, 5) - fy(i, j, k, 5)) / dx[1] +
                         (fz(i, j, k + 1, 5) - fz(i, j, k, 5)) / dx[2]);
                    gp(i, j, k, 1) +=
                        chng * ptfac *
                        ((fx(i + 1, j, k, 6) - fx(i, j, k, 6)) / dx[0] +
                         (fy(i, j + 1, k, 6) - fy(i, j, k, 6)) / dx[1] +
                         (fz(i, j, k + 1, 6) - fz(i, j, k, 6)) / dx[2]);
                    gp(i, j, k, 2) +=
                        chng * ptfac *
                        ((fx(i + 1, j, k, 7) - fx(i, j, k, 7)) / dx[0] +
                         (fy(i, j + 1, k, 7) - fy(i, j, k, 7)) / dx[1] +
                         (fz(i, j, k + 1, 7) - fz(i, j, k, 7)) / dx[2]);
                }

                /*if (gp(i, j, k, 2) > 100.0)
                    std::cout
                        << i << " " << j << " " << k << " gp fluxes "
                        << fx(i + 1, j, k, 7) << " " << -fx(i, j, k, 7) << " "
                        << fy(i, j + 1, k, 7) << " " << -fy(i, j, k, 7) << " "
                        << fz(i, j, k + 1, 7) << " " << -fz(i, j, k, 7)
                        << " density fluxes " << fx(i + 1, j, k, 1) << " "
                        << -fx(i, j, k, 1) << " " << fy(i, j + 1, k, 1) << " "
                        << -fy(i, j, k, 1) << " " << fz(i, j, k + 1, 1) << " "
                        << -fz(i, j, k, 1) << std::endl;*/

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
                const amrex::Real chng = (iblnod(i, j, k) == -1) ? 1.0 : 0.0;
                p(i, j, k) += ptfac * chng * sp(i, j, k);
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

amrex::Real calculate_pseudo_velocity_scale(
    amrex::iMultiFab& mf_iblank,
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dx,
    const amrex::Real pvmax)
{
    // Get minimum length scale from dx
    const amrex::Real dx_min = std::min(std::min(dx[0], dx[1]), dx[2]);
    // The dx of this level should be considered if iblank = -1 here
    // Otherwise, set to max possible value for this mesh (level 0 values)
    const amrex::Real pvscale = (mf_iblank.min(0) == -1) ? dx_min : pvmax;
    return pvscale;
}

amrex::Real calculate_pseudo_dt_flux(
    amrex::MultiFab& mf_fx,
    amrex::MultiFab& mf_fy,
    amrex::MultiFab& mf_fz,
    amrex::MultiFab& mf_vof,
    amrex::iMultiFab& mf_iblank,
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dx,
    amrex::Real tol)
{
    // Get the maximum flux magnitude, but just for vof fluxes
    const amrex::Real pdt_fx = amrex::ReduceMin(
        mf_fx, mf_vof, 0,
        [=] AMREX_GPU_HOST_DEVICE(
            amrex::Box const& bx, amrex::Array4<amrex::Real const> const& fx,
            amrex::Array4<amrex::Real const> const& vof) -> amrex::Real {
            amrex::Real pdt_fab = 1.0;
            amrex::Loop(bx, [=, &pdt_fab](int i, int j, int k) noexcept {
                amrex::Real pdt_lim = 1.0;
                if (fx(i, j, k, 0) > tol && vof(i, j, k) > tol) {
                    // VOF is removed from cell i
                    pdt_lim = vof(i, j, k) * dx[0] / fx(i, j, k, 0);
                } else if (fx(i, j, k, 0) < -tol && vof(i - 1, j, k) > tol) {
                    // VOF is removed from cell i-1
                    pdt_lim = vof(i - 1, j, k) * dx[0] / -fx(i, j, k, 0);
                }
                pdt_fab = amrex::min(pdt_fab, pdt_lim);
            });
            return pdt_fab;
        });
    const amrex::Real pdt_fy = amrex::ReduceMin(
        mf_fy, mf_vof, 0,
        [=] AMREX_GPU_HOST_DEVICE(
            amrex::Box const& bx, amrex::Array4<amrex::Real const> const& fy,
            amrex::Array4<amrex::Real const> const& vof) -> amrex::Real {
            amrex::Real pdt_fab = 1.0;
            amrex::Loop(bx, [=, &pdt_fab](int i, int j, int k) noexcept {
                amrex::Real pdt_lim = 1.0;
                if (fy(i, j, k, 0) > tol && vof(i, j, k) > tol) {
                    // VOF is removed from cell j
                    pdt_lim = vof(i, j, k) * dx[1] / fy(i, j, k, 0);
                } else if (fy(i, j, k, 0) < -tol && vof(i, j - 1, k) > tol) {
                    // VOF is removed from cell j-1
                    pdt_lim = vof(i, j - 1, k) * dx[1] / -fy(i, j, k, 0);
                }
                pdt_fab = amrex::min(pdt_fab, pdt_lim);
            });
            return pdt_fab;
        });
    // Pre-process z fluxes
    for (amrex::MFIter mfi(mf_fz); mfi.isValid(); ++mfi) {
        const auto& vbx = mfi.validbox();
        const amrex::Array4<const int>& iblank = mf_iblank.const_array(mfi);
        const amrex::Array4<amrex::Real>& fz = mf_fz.array(mfi);
        // Populate approximate signed distance function
        amrex::ParallelFor(
            vbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                // Ignore fluxes on overset boundaries
                if (iblank(i, j, k) * iblank(i, j, k - 1) == -1) {
                    fz(i, j, k, 8) = 0.0;
                }
            });
    }
    const amrex::Real pdt_fz = amrex::ReduceMin(
        mf_fz, mf_vof, 0,
        [=] AMREX_GPU_HOST_DEVICE(
            amrex::Box const& bx, amrex::Array4<amrex::Real const> const& fz,
            amrex::Array4<amrex::Real const> const& vof) -> amrex::Real {
            amrex::Real pdt_fab = 1.0;
            amrex::Loop(bx, [=, &pdt_fab](int i, int j, int k) noexcept {
                amrex::Real pdt_lim = 1.0;
                if (fz(i, j, k, 0) * fz(i, j, k, 8) > tol &&
                    vof(i, j, k) > tol) {
                    // VOF is removed from cell k
                    pdt_lim = vof(i, j, k) * dx[2] / fz(i, j, k, 0);
                } else if (
                    fz(i, j, k, 0) * fz(i, j, k, 8) < -tol &&
                    vof(i, j, k - 1) > tol) {
                    // VOF is removed from cell k-1
                    pdt_lim = vof(i, j, k - 1) * dx[2] / -fz(i, j, k, 0);
                }
                pdt_fab = amrex::min(pdt_fab, pdt_lim);
            });
            return pdt_fab;
        });
    const amrex::Real pdt = amrex::min(pdt_fx, amrex::min(pdt_fy, pdt_fz));
    return pdt;
}

// for debugging
void equate_field(amrex::MultiFab& mf_dest, amrex::MultiFab& mf_src)
{
    for (amrex::MFIter mfi(mf_dest); mfi.isValid(); ++mfi) {
        const auto& vbx = mfi.validbox();
        const amrex::Array4<const amrex::Real>& src = mf_src.const_array(mfi);
        const amrex::Array4<amrex::Real>& dest = mf_dest.array(mfi);
        amrex::ParallelFor(
            vbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                dest(i, j, k) = src(i, j, k);
                /*dest(i, j, k) = std::sqrt(
                    src(i, j, k, 0) * src(i, j, k, 0) +
                    src(i, j, k, 1) * src(i, j, k, 1) +
                    src(i, j, k, 2) * src(i, j, k, 2));*/
            });
    }
}

void harmonize_vof(
    amrex::MultiFab& mf_vof_target,
    amrex::MultiFab& mf_vof_original,
    amrex::iMultiFab& mf_iblank)
{
    for (amrex::MFIter mfi(mf_vof_target); mfi.isValid(); ++mfi) {
        const auto& gbx = mfi.growntilebox();
        const amrex::Array4<amrex::Real>& tg_vof = mf_vof_target.array(mfi);
        const amrex::Array4<const amrex::Real>& og_vof =
            mf_vof_original.const_array(mfi);
        const amrex::Array4<const int>& iblank = mf_iblank.const_array(mfi);
        amrex::ParallelFor(
            gbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
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
    amrex::MultiFab& mf_refdens,
    amrex::iMultiFab& mf_iblank,
    const amrex::Real grav_z,
    const bool is_pptb)
{
    for (amrex::MFIter mfi(mf_gp); mfi.isValid(); ++mfi) {
        const auto& gbx = mfi.growntilebox();
        const amrex::Array4<amrex::Real>& gp = mf_gp.array(mfi);
        const amrex::Array4<const amrex::Real>& rho =
            mf_density.const_array(mfi);
        const amrex::Array4<const amrex::Real>& rho0 =
            mf_refdens.const_array(mfi);
        const amrex::Array4<const int>& iblank = mf_iblank.const_array(mfi);
        amrex::ParallelFor(
            gbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                // Replace amr-wind vof values with originals
                if (iblank(i, j, k) == -1) {
                    const amrex::Real dfac =
                        is_pptb ? rho(i, j, k) - rho0(i, j, k) : rho(i, j, k);
                    gp(i, j, k, 0) = 0.;
                    gp(i, j, k, 1) = 0.;
                    gp(i, j, k, 2) = dfac * grav_z;
                }
            });
    }
}

void replace_gradp(
    amrex::MultiFab& mf_gp,
    amrex::MultiFab& mf_gp0,
    amrex::iMultiFab& mf_iblank,
    int ib_val)
{
    for (amrex::MFIter mfi(mf_gp); mfi.isValid(); ++mfi) {
        const auto& gbx = mfi.growntilebox();
        const amrex::Array4<amrex::Real>& gp = mf_gp.array(mfi);
        const amrex::Array4<const amrex::Real>& gp0 = mf_gp0.const_array(mfi);
        const amrex::Array4<const int>& iblank = mf_iblank.const_array(mfi);
        amrex::ParallelFor(
            gbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                // Replace amr-wind vof values with originals
                if (iblank(i, j, k) == ib_val) {
                    gp(i, j, k, 0) = gp0(i, j, k, 0);
                    gp(i, j, k, 1) = gp0(i, j, k, 1);
                    gp(i, j, k, 2) = gp0(i, j, k, 2);
                }
            });
    }
}

void apply_pressure_gradient(
    amrex::MultiFab& mf_vel,
    amrex::MultiFab& mf_density,
    amrex::MultiFab& mf_gp,
    const amrex::Real scaling_factor)
{
    for (amrex::MFIter mfi(mf_gp); mfi.isValid(); ++mfi) {
        const auto& vbx = mfi.validbox();
        const amrex::Array4<amrex::Real>& vel = mf_vel.array(mfi);
        const amrex::Array4<const amrex::Real>& rho =
            mf_density.const_array(mfi);
        const amrex::Array4<const amrex::Real>& gp = mf_gp.const_array(mfi);
        amrex::ParallelFor(
            vbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                const amrex::Real soverrho = scaling_factor / rho(i, j, k);
                vel(i, j, k, 0) -= gp(i, j, k, 0) * soverrho;
                vel(i, j, k, 1) -= gp(i, j, k, 1) * soverrho;
                vel(i, j, k, 2) -= gp(i, j, k, 2) * soverrho;
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
    if (!yesVOF) {
        return;
    }

    // Get phase densities
    auto& mphase = sim.physics_manager().get<MultiPhase>();
    const amrex::Real rhol = mphase.rho1();
    const amrex::Real rhog = mphase.rho2();

    // Get blanking for cells
    auto& iblank_cell = repo.get_int_field("iblank_cell");
    auto& iblank_node = repo.get_int_field("iblank_node");

    // Get fields that will be modified
    auto& vof = repo.get_field("vof");
    auto& levelset = repo.get_field("levelset");
    auto& rho = repo.get_field("density");
    auto& velocity = repo.get_field("velocity");
    auto& gp_noghost = repo.get_field("gp");
    auto& p = repo.get_field("p");

    // Create scratch fields for fluxes - components are vof, density, and
    // 3 of velocity, 3 of pressure gradient, 1 addl component for p_src flag
    auto flux_x = repo.create_scratch_field(9, 1, amr_wind::FieldLoc::XFACE);
    auto flux_y = repo.create_scratch_field(9, 1, amr_wind::FieldLoc::YFACE);
    auto flux_z = repo.create_scratch_field(9, 1, amr_wind::FieldLoc::ZFACE);
    // Create scratch field for pressure source term
    auto p_src = repo.create_scratch_field(1, 0, amr_wind::FieldLoc::NODE);
    // Create scratch fields for normal vector and target vof field
    auto normal_vec = repo.create_scratch_field(3, vof.num_grow()[0] - 1);
    auto target_vof = repo.create_scratch_field(1, vof.num_grow()[0]);
    // Create scratch field for pressure gradient that has ghost cells
    // Sharpening fluxes (at faces) have 1 ghost, requiring fields to have >= 2
    auto gp_scr = repo.create_scratch_field(3, 2);
    auto& gp = *gp_scr;

    // Give initial max possible value of pseudo-velocity scale
    auto dx_lev0 = (geom[0]).CellSizeArray();
    const amrex::Real max_pvscale =
        std::min(std::min(dx_lev0[0], dx_lev0[1]), dx_lev0[2]);
    amrex::Real pvscale = max_pvscale;

    // Prep things that do not change with iterations
    for (int lev = 0; lev < nlevels; ++lev) {
        // Thickness used here is user parameter, whatever works best
        auto dx = (geom[lev]).CellSizeArray();
        const amrex::Real i_th = rlscale * std::cbrt(dx[0] * dx[1] * dx[2]);

        // Populate approximate signed distance function
        populate_psi(levelset(lev), vof(lev), i_th);

        // Populate gp scratch field
        gp(lev).setVal(0.0); // for external boundaries
        amrex::MultiFab::Copy(gp(lev), gp_noghost(lev), 0, 0, 3, 0); // nonghost
        gp(lev).FillBoundary(sim.mesh().Geom()[lev].periodicity());  // internal

        // Get pseudo-velocity scale, proportional to smallest dx in iblank
        const amrex::Real pvscale_lev =
            calculate_pseudo_velocity_scale(iblank_cell(lev), dx, max_pvscale);
        pvscale = std::min(pvscale, pvscale_lev);
    }
    amrex::ParallelDescriptor::ReduceRealMin(pvscale);

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

    // Put fluxes in vector for averaging down during iterations
    amrex::Vector<amrex::Array<amrex::MultiFab*, AMREX_SPACEDIM>> fluxes(
        repo.num_active_levels());
    for (int lev = 0; lev < nlevels; ++lev) {
        fluxes[lev][0] = &(*flux_x)(lev);
        fluxes[lev][1] = &(*flux_y)(lev);
        fluxes[lev][2] = &(*flux_z)(lev);
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

        // Maximum possible value of pseudo time factor (dtau)
        amrex::Real ptfac = 1.0;

        // Maximum pseudoCFL, 0.5 seems to work well
        const amrex::Real pCFL = 0.5;

        for (int lev = 0; lev < nlevels; ++lev) {
            // Populate normal vector
            populate_normal_vector(
                (*normal_vec)(lev), vof(lev), iblank_cell(lev));

            // Sharpening fluxes for vof, density, and momentum
            populate_sharpen_discrete_fluxes(
                iblank_cell(lev), (*flux_x)(lev), (*flux_y)(lev),
                (*flux_z)(lev), vof(lev), (*target_vof)(lev),
                (*normal_vec)(lev), velocity(lev), gp(lev), rho(lev), pvscale,
                margin, rhol, rhog);

            // Process fluxes and get pressure source term
            process_fluxes_calc_src(
                (*flux_x)(lev), (*flux_y)(lev), (*flux_z)(lev), vof(lev),
                (*p_src)(lev), iblank_cell(lev), sh_pressure);

            // Measure convergence to determine if loop can stop
            if (cconv) {
                // Update error at specified interval of steps
                const amrex::Real err_lev =
                    measure_convergence(
                        (*flux_x)(lev), (*flux_y)(lev), (*flux_z)(lev)) /
                    pvscale;
                err = amrex::max(err, err_lev);
            }
        }

        // Average down fluxes across levels for consistency
        for (int lev = nlevels - 1; lev > 0; --lev) {
            amrex::IntVect rr =
                geom[lev].Domain().size() / geom[lev - 1].Domain().size();
            amrex::average_down_faces(
                GetArrOfConstPtrs(fluxes[lev]), fluxes[lev - 1], rr,
                geom[lev - 1]);
        }

        // Get pseudo dt (dtau)
        for (int lev = 0; lev < nlevels; ++lev) {
            // Compare vof fluxes to vof in source cells
            auto dx = (geom[lev]).CellSizeArray();
            // Convergence tolerance determines what size of fluxes matter
            const amrex::Real ptfac_lev = calculate_pseudo_dt_flux(
                (*flux_x)(lev), (*flux_y)(lev), (*flux_z)(lev), vof(lev),
                iblank_cell(lev), dx, tol);
            ptfac = amrex::min(ptfac, ptfac_lev);
        }
        amrex::ParallelDescriptor::ReduceRealMin(ptfac);

        // Conform pseudo dt (dtau) to pseudo CFL
        ptfac = pCFL * ptfac;

        // Apply fluxes
        for (int lev = 0; lev < nlevels; ++lev) {
            auto dx = (geom[lev]).CellSizeArray();

            apply_fluxes(
                iblank_cell(lev), iblank_node(lev), (*flux_x)(lev),
                (*flux_y)(lev), (*flux_z)(lev), (*p_src)(lev), vof(lev),
                rho(lev), velocity(lev), gp(lev), p(lev), dx, ptfac,
                sh_pressure);

            // Update ghost cells
            vof(lev).FillBoundary(sim.mesh().Geom()[lev].periodicity());
            velocity(lev).FillBoundary(sim.mesh().Geom()[lev].periodicity());
            gp(lev).FillBoundary(sim.mesh().Geom()[lev].periodicity());
        }

        // Update density according to vof (fillpatch built in)
        mphase.set_density_via_vof();

        // Copy gp so it can be seen
        amr_wind::field_ops::copy(
            gp_noghost, gp, 0, 0, gp.num_comp(), gp_noghost.num_grow());
        // Output for debugging
        // sim.io_manager().write_plot_file();
        // sim.time().increment_timestep();

        // Ensure that err is same across processors
        if (cconv) {
            amrex::ParallelDescriptor::ReduceRealMax(err);
        }

        amrex::Print() << "sharpen step " << n << " " << err << " " << tol
                       << std::endl;
    }

    // Purely for debugging via visualization, should be removed later
    // Currently set up to overwrite the levelset field (not used as time
    // evolves) with the post-sharpening velocity magnitude
    for (int lev = 0; lev < nlevels; ++lev) {
        equate_field(levelset(lev), (*target_vof)(lev));
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
    const bool is_pptb = mphase.perturb_p();

    // Get blanking for cells
    auto& iblank_cell = repo.get_int_field("iblank_cell");

    // Get fields that will be modified or used
    Field* rho0{nullptr};
    auto& rho = repo.get_field("density");
    auto& gp = repo.get_field("gp");
    if (is_pptb) {
        rho0 = &(sim.repo().get_field("reference_density"));
    } else {
        // Point to existing field, won't be used
        rho0 = &rho;
    }

    // Replace initial gp with best guess (hydrostatic)
    for (int lev = 0; lev < nlevels; ++lev) {
        replace_gradp_hs(
            gp(lev), rho(lev), (*rho0)(lev), iblank_cell(lev), grav_z, is_pptb);
    }
}

void ReapplyModifiedGradP(CFDSim& sim)
{
    auto& repo = sim.repo();
    auto nlevels = repo.num_active_levels();
    auto geom = sim.mesh().Geom();

    // Get gravity
    auto& mphase = sim.physics_manager().get<MultiPhase>();
    const amrex::Real grav_z = mphase.gravity()[2];
    const bool is_pptb = mphase.perturb_p();

    // Get timestep
    const amrex::Real dt = sim.time().deltaT();

    // Get blanking for cells
    auto& iblank_cell = repo.get_int_field("iblank_cell");

    // Get fields that will be modified or used
    Field* rho0{nullptr};
    auto& vel = repo.get_field("velocity");
    auto& rho = repo.get_field("density");
    auto& gp = repo.get_field("gp");
    if (is_pptb) {
        rho0 = &(sim.repo().get_field("reference_density"));
    } else {
        // Point to existing field, won't be used
        rho0 = &rho;
    }

    // For iblanked cells, replace gp with original gp, to get original vel
    for (int lev = 0; lev < nlevels; ++lev) {
        // Remove pressure gradient term
        apply_pressure_gradient(vel(lev), rho(lev), gp(lev), -dt);
        // Modify pressure gradient
        replace_gradp_hs(
            gp(lev), rho(lev), (*rho0)(lev), iblank_cell(lev), grav_z, is_pptb);
        // Reapply pressure gradient term
        apply_pressure_gradient(vel(lev), rho(lev), gp(lev), dt);
    }
}

void CopyGradP(ScratchField& gp_copy, Field& gp, const int nlevels)
{
    // For iblanked cells, replace gp with original gp, to get original vel
    for (int lev = 0; lev < nlevels; ++lev) {
        amrex::MultiFab::Copy(
            gp_copy(lev), gp(lev), 0, 0, gp(lev).nComp(), gp(lev).nGrow());
    }
}

void ReapplyOversetGradP(ScratchField& gp_copy, CFDSim& sim)
{
    auto& repo = sim.repo();
    auto nlevels = repo.num_active_levels();

    // Get timestep
    const amrex::Real dt = sim.time().deltaT();

    // Get blanking for cells
    auto& iblank_cell = repo.get_int_field("iblank_cell");

    // Get fields that will be modified or used
    auto& vel = repo.get_field("velocity");
    auto& rho = repo.get_field("density");
    auto& gp = repo.get_field("gp");

    // For iblanked cells, replace gp with original gp, to get original vel
    for (int lev = 0; lev < nlevels; ++lev) {
        // Remove pressure gradient term
        apply_pressure_gradient(vel(lev), rho(lev), gp(lev), -dt);
        // Modify pressure gradient, replace blanked version
        replace_gradp(gp(lev), gp_copy(lev), iblank_cell(lev), -1);
        // Reapply pressure gradient term
        apply_pressure_gradient(vel(lev), rho(lev), gp(lev), dt);
    }
}

void ReplaceUnMaskedGradP(
    Field& gp, ScratchField& gp_copy, IntField& iblank_cell)
{
    // Replace gp with copy from before exchange
    for (int lev = 0; lev < gp.repo().num_active_levels(); ++lev) {
        replace_gradp(gp(lev), gp_copy(lev), iblank_cell(lev), 1);
    }
}

} // namespace amr_wind::overset

#endif