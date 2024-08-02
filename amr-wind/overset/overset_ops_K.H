#ifndef OVERSET_OPS_K_H_
#define OVERSET_OPS_K_H_

#include "amr-wind/core/vs/vector_space.H"
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
    const amrex::IntVect dv{(int)(dir == 0), (int)(dir == 1), (int)(dir == 2)};
    const amrex::IntVect ivm = iv - dv;

    // Gradient of phi normal to interface
    const amrex::Real gphi = (vof(iv) - vof(ivm));
    // Normal vector in each cell (already normalized)
    const amrex::Real norm_ = std::abs(normal(iv, dir));
    const amrex::Real norm_nb = std::abs(normal(ivm, dir));

    // Determine which delta_phi (and multiply by normal)
    // The sign depends on side of flux face
    const amrex::Real dphi_ = (vof(iv) - tg_vof(iv)) * norm_;
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
    const amrex::IntVect dv{(int)(dir == 0), (int)(dir == 1), (int)(dir == 2)};
    const amrex::IntVect ivm = iv - dv;

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

void AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE gp_rho_face(
    const int i,
    const int j,
    const int k,
    const int dir,
    amrex::Array4<amrex::Real const> const& vof,
    amrex::Array4<amrex::Real const> const& gp,
    amrex::Array4<amrex::Real const> const& rho,
    amrex::Real& uface,
    amrex::Real& vface,
    amrex::Real& wface)
{
    // Set up neighbor indices
    const amrex::IntVect iv{i, j, k};
    const amrex::IntVect dv{(int)(dir == 0), (int)(dir == 1), (int)(dir == 2)};
    const amrex::IntVect ivm = iv - dv;

    // Gradient of phi normal to interface
    const amrex::Real gphi = (vof(iv) - vof(ivm));

    // Get velocities on both sides
    const amrex::Real u_ = gp(iv, 0) / rho(iv);
    const amrex::Real v_ = gp(iv, 1) / rho(iv);
    const amrex::Real w_ = gp(iv, 2) / rho(iv);
    const amrex::Real u_nb = gp(ivm, 0) / rho(ivm);
    const amrex::Real v_nb = gp(ivm, 1) / rho(ivm);
    const amrex::Real w_nb = gp(ivm, 2) / rho(ivm);
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

vs::Tensor AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE gp_flux_tensor(
    const int i,
    const int j,
    const int k,
    amrex::Array4<amrex::Real const> const& fx,
    amrex::Array4<amrex::Real const> const& fy,
    amrex::Array4<amrex::Real const> const& fz,
    const amrex::Real tiny)
{
    // Averaging depends on iblank array, encoded in 8th component
    const amrex::Real avg_fx = fx(i, j, k, 8) + fx(i, j, k - 1, 8) +
                               fx(i, j - 1, k, 8) + fx(i, j - 1, k - 1, 8) +
                               tiny;
    const amrex::Real avg_fy = fy(i, j, k, 8) + fy(i - 1, j, k, 8) +
                               fy(i, j, k - 1, 8) + fy(i - 1, j, k - 1, 8) +
                               tiny;
    const amrex::Real avg_fz = fz(i, j, k, 8) + fz(i - 1, j, k, 8) +
                               fz(i, j - 1, k, 8) + fz(i - 1, j - 1, k, 8) +
                               tiny;

    auto f_otimes_gradp = vs::Tensor::identity();
    // Average fluxes (faces) to nodes where pressure exists
    f_otimes_gradp.xx() = (fx(i, j, k, 5) + fx(i, j, k - 1, 5) +
                           fx(i, j - 1, k, 5) + fx(i, j - 1, k - 1, 5)) /
                          avg_fx;
    f_otimes_gradp.xy() = (fx(i, j, k, 6) + fx(i, j, k - 1, 6) +
                           fx(i, j - 1, k, 6) + fx(i, j - 1, k - 1, 6)) /
                          avg_fx;
    f_otimes_gradp.xz() = (fx(i, j, k, 7) + fx(i, j, k - 1, 7) +
                           fx(i, j - 1, k, 7) + fx(i, j - 1, k - 1, 7)) /
                          avg_fx;
    f_otimes_gradp.yx() = (fy(i, j, k, 5) + fy(i - 1, j, k, 5) +
                           fy(i, j, k - 1, 5) + fy(i - 1, j, k - 1, 5)) /
                          avg_fy;
    f_otimes_gradp.yy() = (fy(i, j, k, 6) + fy(i - 1, j, k, 6) +
                           fy(i, j, k - 1, 6) + fy(i - 1, j, k - 1, 6)) /
                          avg_fy;
    f_otimes_gradp.yz() = (fy(i, j, k, 7) + fy(i - 1, j, k, 7) +
                           fy(i, j, k - 1, 7) + fy(i - 1, j, k - 1, 7)) /
                          avg_fy;
    f_otimes_gradp.zx() = (fz(i, j, k, 5) + fz(i - 1, j, k, 5) +
                           fz(i, j - 1, k, 5) + fz(i - 1, j - 1, k, 5)) /
                          avg_fz;
    f_otimes_gradp.zy() = (fz(i, j, k, 6) + fz(i - 1, j, k, 6) +
                           fz(i, j - 1, k, 6) + fz(i - 1, j - 1, k, 6)) /
                          avg_fz;
    f_otimes_gradp.zz() = (fz(i, j, k, 7) + fz(i - 1, j, k, 7) +
                           fz(i, j - 1, k, 7) + fz(i - 1, j - 1, k, 7)) /
                          avg_fz;

    return f_otimes_gradp;
}

vs::Tensor AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE normal_reinit_tensor(
    const int i,
    const int j,
    const int k,
    amrex::Array4<amrex::Real const> const& fx,
    amrex::Array4<amrex::Real const> const& fy,
    amrex::Array4<amrex::Real const> const& fz,
    amrex::Array4<amrex::Real const> const& vof,
    const amrex::Real tiny)
{
    // Averaging depends on iblank array, encoded in 8th component
    const amrex::Real avg_fx = fx(i, j, k, 8) + fx(i, j, k - 1, 8) +
                               fx(i, j - 1, k, 8) + fx(i, j - 1, k - 1, 8) +
                               tiny;
    const amrex::Real avg_fy = fy(i, j, k, 8) + fy(i - 1, j, k, 8) +
                               fy(i, j, k - 1, 8) + fy(i - 1, j, k - 1, 8) +
                               tiny;
    const amrex::Real avg_fz = fz(i, j, k, 8) + fz(i - 1, j, k, 8) +
                               fz(i, j - 1, k, 8) + fz(i - 1, j - 1, k, 8) +
                               tiny;

    auto n_zeta = vs::Vector::one();
    // Get components of normal in each position
    n_zeta.x() =
        ((vof(i, j, k) - vof(i - 1, j, k)) * fx(i, j, k, 8) +
         (vof(i, j - 1, k) - vof(i - 1, j - 1, k)) * fx(i, j - 1, k, 8) +
         (vof(i, j, k - 1) - vof(i - 1, j, k - 1)) * fx(i, j, k - 1, 8) +
         (vof(i, j - 1, k - 1) - vof(i - 1, j - 1, k - 1)) *
             fx(i, j - 1, k - 1, 8)) /
        avg_fx;
    n_zeta.y() =
        ((vof(i, j, k) - vof(i, j - 1, k)) * fy(i, j, k, 8) +
         (vof(i - 1, j, k) - vof(i - 1, j - 1, k)) * fy(i - 1, j, k, 8) +
         (vof(i, j, k - 1) - vof(i, j - 1, k - 1)) * fy(i, j, k - 1, 8) +
         (vof(i - 1, j, k - 1) - vof(i - 1, j - 1, k - 1)) *
             fy(i - 1, j, k - 1, 8)) /
        avg_fy;
    n_zeta.z() =
        ((vof(i, j, k) - vof(i, j, k - 1)) * fz(i, j, k, 8) +
         (vof(i - 1, j, k) - vof(i - 1, j, k - 1)) * fz(i - 1, j, k, 8) +
         (vof(i, j - 1, k) - vof(i, j - 1, k - 1)) * fz(i, j - 1, k, 8) +
         (vof(i - 1, j - 1, k) - vof(i - 1, j - 1, k - 1)) *
             fz(i - 1, j - 1, k, 8)) /
        avg_fz;

    n_zeta.normalize();

    // To do outer product between vectors, set up tensors
    auto n_tensor = vs::Tensor::identity();
    auto n_tensor_T = vs::Tensor::identity();
    n_tensor.rows(n_zeta, n_zeta, n_zeta);
    n_tensor_T.cols(n_zeta, n_zeta, n_zeta);
    // Multiply tensors elementwise
    auto n_otimes_n = vs::Tensor::identity();
    n_otimes_n.rows(
        n_tensor.x() * n_tensor_T.x(), n_tensor.y() * n_tensor_T.y(),
        n_tensor.z() * n_tensor_T.z());

    return n_otimes_n;
}

} // namespace amr_wind::overset_ops

#endif