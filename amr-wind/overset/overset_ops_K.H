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
    const int i,
    const int j,
    const int k,
    const int dir,
    const amrex::Real margin,
    amrex::Array4<amrex::Real const> const& vof,
    amrex::Array4<amrex::Real const> const& tg_vof,
    amrex::Array4<amrex::Real const> const& normal)
{
    // Set up neighbor indices
    const amrex::IntVect iv{i, j, k};
    const amrex::IntVect ivm =
        iv -
        (amrex::IntVect){(int)(dir == 0), (int)(dir == 1), (int)(dir == 2)};

    // Gradient of phi normal to interface
    const amrex::Real gphi = (vof(iv) - vof(ivm));
    // Normal vector in each cell (already normalized)
    const amrex::Real norm_ = normal(iv, dir);
    const amrex::Real norm_nb = normal(ivm, dir);

    // Determine which delta_phi (and multiply by normal)
    // The sign depends on side of flux face (like upwinding)
    const amrex::Real dphi_ = (tg_vof(iv) - vof(iv)) * (-norm_);
    const amrex::Real dphi_nb = (tg_vof(ivm) - vof(ivm)) * norm_nb;
    // Average value used across the interface
    amrex::Real dphi_eval = 0.5 * (dphi_ + dphi_nb);
    // Upwinding when on the gas side, downwinding on the liquid
    // Across the interface defined as crossing 0.5 or within margin of 0.5
    if ((std::abs(vof(iv) - 0.5) > margin ||
         std::abs(vof(ivm) - 0.5) > margin)) {
        if (gphi > 0.0) {
            dphi_eval = (vof(ivm) < 0.5 && vof(iv) <= 0.5 + margin) ? dphi_nb
                                                                    : dphi_eval;
            dphi_eval =
                (vof(ivm) >= 0.5 - margin && vof(iv) > 0.5) ? dphi_ : dphi_eval;
        }
        if (gphi < 0.0) {
            dphi_eval =
                (vof(iv) < 0.5 && vof(ivm) <= 0.5 + margin) ? dphi_ : dphi_eval;
            dphi_eval = (vof(iv) >= 0.5 - margin && vof(ivm) > 0.5) ? dphi_nb
                                                                    : dphi_eval;
        }
    }
    return dphi_eval;
}

void AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE velocity_face(
    const int i,
    const int j,
    const int k,
    const int dir,
    amrex::Array4<amrex::Real const> const& vof,
    amrex::Array4<amrex::Real const> const& velocity,
    amrex::Real& uface,
    amrex::Real& vface,
    amrex::Real& wface)
{
    // Set up neighbor indices
    const amrex::IntVect iv{i, j, k};
    const amrex::IntVect ivm =
        iv -
        (amrex::IntVect){(int)(dir == 0), (int)(dir == 1), (int)(dir == 2)};

    // Gradient of phi normal to interface
    const amrex::Real gphi = (vof(iv) - vof(ivm));

    // Get velocities on both sides
    const amrex::Real u_ = velocity(iv, 0);
    const amrex::Real v_ = velocity(iv, 1);
    const amrex::Real w_ = velocity(iv, 2);
    const amrex::Real u_nb = velocity(ivm, 0);
    const amrex::Real v_nb = velocity(ivm, 1);
    const amrex::Real w_nb = velocity(ivm, 2);
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