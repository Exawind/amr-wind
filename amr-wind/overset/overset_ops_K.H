#ifndef OVERSET_OPS_K_H_
#define OVERSET_OPS_K_H_

#include <AMReX_FArrayBox.H>

namespace amr_wind::overset_ops {

// Approximate signed distance function
amrex::Real AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
asdf(const amrex::Real a_vof, const amrex::Real i_th, const amrex::Real tiny)
{
    // function of local vof value and interface thickness
    return (i_th * log((a_vof + tiny) / (1. - a_vof + tiny)));
}

amrex::Real AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE alpha_flux(
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

void AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE velocity_face(
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

} // namespace amr_wind::overset_ops

#endif