#ifndef OVERSET_OPS_ROUTINES_H_
#define OVERSET_OPS_ROUTINES_H_

#include "amr-wind/equation_systems/vof/volume_fractions.H"
#include "amr-wind/overset/overset_ops_K.H"

namespace amr_wind::overset_ops {

// Populate approximate signed distance function using vof field
void populate_psi(
    amrex::MultiFab& mf_psi,
    amrex::MultiFab& mf_vof,
    const amrex::Real i_th,
    const amrex::Real asdf_tiny)
{
    for (amrex::MFIter mfi(mf_psi); mfi.isValid(); ++mfi) {
        const auto& gbx = mfi.growntilebox();
        const amrex::Array4<amrex::Real>& psi = mf_psi.array(mfi);
        const amrex::Array4<const amrex::Real>& vof = mf_vof.const_array(mfi);
        amrex::ParallelFor(
            gbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                psi(i, j, k) = asdf(vof(i, j, k), i_th, asdf_tiny);
            });
    }
}

// Modify a vof field to not have values that barely differ from 0 or 1
void process_vof(amrex::MultiFab& mf_vof, const amrex::Real vof_tol)
{
    for (amrex::MFIter mfi(mf_vof); mfi.isValid(); ++mfi) {
        const auto& gbx = mfi.growntilebox();
        const amrex::Array4<amrex::Real>& vof = mf_vof.array(mfi);
        amrex::ParallelFor(
            gbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                vof(i, j, k) =
                    vof(i, j, k) < vof_tol
                        ? 0.0
                        : (vof(i, j, k) > 1. - vof_tol ? 1. : vof(i, j, k));
            });
    }
}

// Combine overset target vof field with current non-overset vof field
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

// Populate normal vector with special treatment of overset boundary
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

// Calculate fluxes for reinitialization over entire domain without concern for
// overset bdy
void populate_sharpen_fluxes(
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
        const amrex::Array4<amrex::Real>& fx = mf_fx.array(mfi);
        const amrex::Array4<amrex::Real>& fy = mf_fy.array(mfi);
        const amrex::Array4<amrex::Real>& fz = mf_fz.array(mfi);
        const amrex::Array4<const amrex::Real>& vof = mf_vof.const_array(mfi);
        const amrex::Array4<const amrex::Real>& tg_vof =
            mf_target_vof.const_array(mfi);
        const amrex::Array4<const amrex::Real>& norm = mf_norm.const_array(mfi);
        const amrex::Array4<const amrex::Real>& vel =
            mf_velocity.const_array(mfi);
        const amrex::Array4<const amrex::Real>& gp = mf_gp.const_array(mfi);
        const amrex::Array4<const amrex::Real>& rho =
            mf_density.const_array(mfi);
        // Populate vof and density fluxes for each direction
        amrex::ParallelFor(
            xbxp1, ybxp1, zbxp1,
            [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                // vof flux
                amrex::Real flux =
                    Gamma * alpha_flux(i, j, k, 0, margin, vof, tg_vof, norm);
                fx(i, j, k, 0) = flux;
                // density flux
                flux *= (rho1 - rho2);
                fx(i, j, k, 1) = flux;
                // momentum fluxes (dens flux * face vel)
                amrex::Real uf, vf, wf;
                velocity_face(i, j, k, 0, vof, vel, uf, vf, wf);
                fx(i, j, k, 2) = flux * uf;
                fx(i, j, k, 3) = flux * vf;
                fx(i, j, k, 4) = flux * wf;
                gp_rho_face(i, j, k, 0, vof, gp, rho, uf, vf, wf);
                fx(i, j, k, 5) = flux * uf;
                fx(i, j, k, 6) = flux * vf;
                fx(i, j, k, 7) = flux * wf;
                // Turn "on" all flux faces, later modified in
                // process_fluxes_calc_src
                fx(i, j, k, 8) = 1.0;
            },
            [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                amrex::Real flux =
                    Gamma * alpha_flux(i, j, k, 1, margin, vof, tg_vof, norm);
                fy(i, j, k, 0) = flux;
                flux *= (rho1 - rho2);
                fy(i, j, k, 1) = flux;
                amrex::Real uf, vf, wf;
                velocity_face(i, j, k, 1, vof, vel, uf, vf, wf);
                fy(i, j, k, 2) = flux * uf;
                fy(i, j, k, 3) = flux * vf;
                fy(i, j, k, 4) = flux * wf;
                gp_rho_face(i, j, k, 1, vof, gp, rho, uf, vf, wf);
                fy(i, j, k, 5) = flux * uf;
                fy(i, j, k, 6) = flux * vf;
                fy(i, j, k, 7) = flux * wf;
                fy(i, j, k, 8) = 1.0;
            },
            [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                amrex::Real flux =
                    Gamma * alpha_flux(i, j, k, 2, margin, vof, tg_vof, norm);
                fz(i, j, k, 0) = flux;
                flux *= (rho1 - rho2);
                fz(i, j, k, 1) = flux;
                amrex::Real uf, vf, wf;
                velocity_face(i, j, k, 2, vof, vel, uf, vf, wf);
                fz(i, j, k, 2) = flux * uf;
                fz(i, j, k, 3) = flux * vf;
                fz(i, j, k, 4) = flux * wf;
                gp_rho_face(i, j, k, 2, vof, gp, rho, uf, vf, wf);
                fz(i, j, k, 5) = flux * uf;
                fz(i, j, k, 6) = flux * vf;
                fz(i, j, k, 7) = flux * wf;
                fz(i, j, k, 8) = 1.0;
            });
    }
}

// Process reinitialization fluxes - zero non-internal to overset region;
// also calculate pressure source / sink term as a function of fluxes
void process_fluxes_calc_src(
    amrex::MultiFab& mf_fx,
    amrex::MultiFab& mf_fy,
    amrex::MultiFab& mf_fz,
    amrex::MultiFab& mf_vof,
    amrex::MultiFab& mf_psource,
    amrex::iMultiFab& mf_iblank)
{
    int ncompx = mf_fx.nComp();
    int ncompy = mf_fy.nComp();
    int ncompz = mf_fz.nComp();
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
                bool zero_all = (iblank(i, j, k - 1) + iblank(i, j, k) > -2);
                for (int n = 0; n < ncompz; ++n) {
                    fz(i, j, k, n) *= zero_all ? 0. : 1.;
                }
            });
        // With knowledge of fluxes, compute pressure source term
        amrex::Box const& nbx = mfi.nodaltilebox();
        constexpr amrex::Real tiny =
            std::numeric_limits<amrex::Real>::epsilon();
        amrex::ParallelFor(
            nbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                // Get components of outer product between f and gradp/rho
                // from interpolating between values at neighboring faces
                const amrex::Real fxgpx_flux =
                    (fx(i, j, k, 5) + fx(i, j, k - 1, 5) + fx(i, j - 1, k, 5) +
                     fx(i, j - 1, k - 1, 5)) /
                    (fx(i, j, k, 8) + fx(i, j, k - 1, 8) + fx(i, j - 1, k, 8) +
                     fx(i, j - 1, k - 1, 8) + tiny);
                const amrex::Real fxgpy_flux =
                    (fx(i, j, k, 6) + fx(i, j, k - 1, 6) + fx(i, j - 1, k, 6) +
                     fx(i, j - 1, k - 1, 6)) /
                    (fx(i, j, k, 8) + fx(i, j, k - 1, 8) + fx(i, j - 1, k, 8) +
                     fx(i, j - 1, k - 1, 8) + tiny);
                const amrex::Real fxgpz_flux =
                    (fx(i, j, k, 7) + fx(i, j, k - 1, 7) + fx(i, j - 1, k, 7) +
                     fx(i, j - 1, k - 1, 7)) /
                    (fx(i, j, k, 8) + fx(i, j, k - 1, 8) + fx(i, j - 1, k, 8) +
                     fx(i, j - 1, k - 1, 8) + tiny);

                const amrex::Real fygpx_flux =
                    (fy(i, j, k, 5) + fy(i - 1, j, k, 5) + fy(i, j, k - 1, 5) +
                     fy(i - 1, j, k - 1, 5)) /
                    (fy(i, j, k, 8) + fy(i - 1, j, k, 8) + fy(i, j, k - 1, 8) +
                     fy(i - 1, j, k - 1, 8) + tiny);
                const amrex::Real fygpy_flux =
                    (fy(i, j, k, 6) + fy(i - 1, j, k, 6) + fy(i, j, k - 1, 6) +
                     fy(i - 1, j, k - 1, 6)) /
                    (fy(i, j, k, 8) + fy(i - 1, j, k, 8) + fy(i, j, k - 1, 8) +
                     fy(i - 1, j, k - 1, 8) + tiny);
                const amrex::Real fygpz_flux =
                    (fy(i, j, k, 7) + fy(i - 1, j, k, 7) + fy(i, j, k - 1, 7) +
                     fy(i - 1, j, k - 1, 7)) /
                    (fy(i, j, k, 8) + fy(i - 1, j, k, 8) + fy(i, j, k - 1, 8) +
                     fy(i - 1, j, k - 1, 8) + tiny);

                const amrex::Real fzgpx_flux =
                    (fz(i, j, k, 5) + fz(i - 1, j, k, 5) + fz(i, j - 1, k, 5) +
                     fz(i - 1, j - 1, k, 5)) /
                    (fz(i, j, k, 8) + fz(i - 1, j, k, 8) + fz(i, j - 1, k, 8) +
                     fz(i - 1, j - 1, k, 8) + tiny);
                const amrex::Real fzgpy_flux =
                    (fz(i, j, k, 6) + fz(i - 1, j, k, 6) + fz(i, j - 1, k, 6) +
                     fz(i - 1, j - 1, k, 6)) /
                    (fz(i, j, k, 8) + fz(i - 1, j, k, 8) + fz(i, j - 1, k, 8) +
                     fz(i - 1, j - 1, k, 8) + tiny);
                const amrex::Real fzgpz_flux =
                    (fz(i, j, k, 7) + fz(i - 1, j, k, 7) + fz(i, j - 1, k, 7) +
                     fz(i - 1, j - 1, k, 7)) /
                    (fz(i, j, k, 8) + fz(i - 1, j, k, 8) + fz(i, j - 1, k, 8) +
                     fz(i - 1, j - 1, k, 8) + tiny);

                // Calculate interface normal at pressure node
                amrex::Real n0 =
                    ((vof(i, j, k) - vof(i - 1, j, k)) * fx(i, j, k, 8) +
                     (vof(i, j - 1, k) - vof(i - 1, j - 1, k)) *
                         fx(i, j - 1, k, 8) +
                     (vof(i, j, k - 1) - vof(i - 1, j, k - 1)) *
                         fx(i, j, k - 1, 8) +
                     (vof(i, j - 1, k - 1) - vof(i - 1, j - 1, k - 1)) *
                         fx(i, j - 1, k - 1, 8)) /
                    (fx(i, j, k, 8) + fx(i, j, k - 1, 8) + fx(i, j - 1, k, 8) +
                     fx(i, j - 1, k - 1, 8) + tiny);
                amrex::Real n1 =
                    ((vof(i, j, k) - vof(i, j - 1, k)) * fy(i, j, k, 8) +
                     (vof(i - 1, j, k) - vof(i - 1, j - 1, k)) *
                         fy(i - 1, j, k, 8) +
                     (vof(i, j, k - 1) - vof(i, j - 1, k - 1)) *
                         fy(i, j, k - 1, 8) +
                     (vof(i - 1, j, k - 1) - vof(i - 1, j - 1, k - 1)) *
                         fy(i - 1, j, k - 1, 8)) /
                    (fy(i, j, k, 8) + fy(i, j, k - 1, 8) + fy(i - 1, j, k, 8) +
                     fy(i - 1, j, k - 1, 8) + tiny);
                amrex::Real n2 =
                    ((vof(i, j, k) - vof(i, j, k - 1)) * fz(i, j, k, 8) +
                     (vof(i - 1, j, k) - vof(i - 1, j, k - 1)) *
                         fz(i - 1, j, k, 8) +
                     (vof(i, j - 1, k) - vof(i, j - 1, k - 1)) *
                         fz(i, j - 1, k, 8) +
                     (vof(i - 1, j - 1, k) - vof(i - 1, j - 1, k - 1)) *
                         fz(i - 1, j - 1, k, 8)) /
                    (fz(i, j, k, 8) + fz(i, j - 1, k, 8) + fz(i - 1, j, k, 8) +
                     fz(i - 1, j - 1, k, 8) + tiny);
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

// Calculate a type of CFL by measuring how much % VOF is being removed per cell
amrex::Real calculate_pseudo_dt_flux(
    amrex::MultiFab& mf_fx,
    amrex::MultiFab& mf_fy,
    amrex::MultiFab& mf_fz,
    amrex::MultiFab& mf_vof,
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
                    pdt_lim = vof(i, j, k) / fx(i, j, k, 0);
                } else if (fx(i, j, k, 0) < -tol && vof(i - 1, j, k) > tol) {
                    // VOF is removed from cell i-1
                    pdt_lim = vof(i - 1, j, k) / -fx(i, j, k, 0);
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
                    pdt_lim = vof(i, j, k) / fy(i, j, k, 0);
                } else if (fy(i, j, k, 0) < -tol && vof(i, j - 1, k) > tol) {
                    // VOF is removed from cell j-1
                    pdt_lim = vof(i, j - 1, k) / -fy(i, j, k, 0);
                }
                pdt_fab = amrex::min(pdt_fab, pdt_lim);
            });
            return pdt_fab;
        });
    const amrex::Real pdt_fz = amrex::ReduceMin(
        mf_fz, mf_vof, 0,
        [=] AMREX_GPU_HOST_DEVICE(
            amrex::Box const& bx, amrex::Array4<amrex::Real const> const& fz,
            amrex::Array4<amrex::Real const> const& vof) -> amrex::Real {
            amrex::Real pdt_fab = 1.0;
            amrex::Loop(bx, [=, &pdt_fab](int i, int j, int k) noexcept {
                amrex::Real pdt_lim = 1.0;
                if (fz(i, j, k, 0) > tol && vof(i, j, k) > tol) {
                    // VOF is removed from cell k
                    pdt_lim = vof(i, j, k) / fz(i, j, k, 0);
                } else if (fz(i, j, k, 0) < -tol && vof(i, j, k - 1) > tol) {
                    // VOF is removed from cell k-1
                    pdt_lim = vof(i, j, k - 1) / -fz(i, j, k, 0);
                }
                pdt_fab = amrex::min(pdt_fab, pdt_lim);
            });
            return pdt_fab;
        });
    const amrex::Real pdt = amrex::min(pdt_fx, amrex::min(pdt_fy, pdt_fz));
    return pdt;
}

// Apply reinitialization fluxes to modify fields
void apply_fluxes(
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
    amrex::Real vof_tol)
{

    for (amrex::MFIter mfi(mf_vof); mfi.isValid(); ++mfi) {
        const auto& vbx = mfi.validbox();
        const auto& nbx = mfi.nodaltilebox();
        const amrex::Array4<amrex::Real>& fx = mf_fx.array(mfi);
        const amrex::Array4<amrex::Real>& fy = mf_fy.array(mfi);
        const amrex::Array4<amrex::Real>& fz = mf_fz.array(mfi);
        const amrex::Array4<amrex::Real>& vof = mf_vof.array(mfi);
        const amrex::Array4<amrex::Real>& dens = mf_dens.array(mfi);
        const amrex::Array4<amrex::Real>& vel = mf_vel.array(mfi);
        const amrex::Array4<amrex::Real>& gp = mf_gp.array(mfi);
        const amrex::Array4<amrex::Real>& p = mf_pressure.array(mfi);
        const amrex::Array4<amrex::Real>& sp = mf_psource.array(mfi);
        amrex::ParallelFor(
            vbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
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

                gp(i, j, k, 0) +=
                    ptfac * ((fx(i + 1, j, k, 5) - fx(i, j, k, 5)) / dx[0] +
                             (fy(i, j + 1, k, 5) - fy(i, j, k, 5)) / dx[1] +
                             (fz(i, j, k + 1, 5) - fz(i, j, k, 5)) / dx[2]);
                gp(i, j, k, 1) +=
                    ptfac * ((fx(i + 1, j, k, 6) - fx(i, j, k, 6)) / dx[0] +
                             (fy(i, j + 1, k, 6) - fy(i, j, k, 6)) / dx[1] +
                             (fz(i, j, k + 1, 6) - fz(i, j, k, 6)) / dx[2]);
                gp(i, j, k, 2) +=
                    ptfac * ((fx(i + 1, j, k, 7) - fx(i, j, k, 7)) / dx[0] +
                             (fy(i, j + 1, k, 7) - fy(i, j, k, 7)) / dx[1] +
                             (fz(i, j, k + 1, 7) - fz(i, j, k, 7)) / dx[2]);

                // Ensure vof is bounded
                vof(i, j, k) =
                    vof(i, j, k) < vof_tol
                        ? 0.0
                        : (vof(i, j, k) > 1. - vof_tol ? 1. : vof(i, j, k));
                // Density bounds are enforced elsewhere
            });

        // Apply pressure source/sink
        amrex::ParallelFor(
            nbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                p(i, j, k) += ptfac * sp(i, j, k);
            });
    }
}

// Get the size of the smallest VOF flux to quantify convergence
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

// Set levelset field to another quantity to view in plotfile for debugging
void equate_field(amrex::MultiFab& mf_dest, amrex::MultiFab& mf_src)
{
    for (amrex::MFIter mfi(mf_dest); mfi.isValid(); ++mfi) {
        const auto& vbx = mfi.validbox();
        const amrex::Array4<const amrex::Real>& src = mf_src.const_array(mfi);
        const amrex::Array4<amrex::Real>& dest = mf_dest.array(mfi);
        amrex::ParallelFor(
            vbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                dest(i, j, k) = std::sqrt(
                    src(i, j, k, 0) * src(i, j, k, 0) +
                    src(i, j, k, 1) * src(i, j, k, 1) +
                    src(i, j, k, 2) * src(i, j, k, 2));
            });
    }
}

// Replace pressure gradient with hydrostatic field in overset regions
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

// Swap pressure gradient values in overset region
void replace_gradp(
    amrex::MultiFab& mf_gp,
    amrex::MultiFab& mf_gp0,
    amrex::iMultiFab& mf_iblank)
{
    for (amrex::MFIter mfi(mf_gp); mfi.isValid(); ++mfi) {
        const auto& gbx = mfi.growntilebox();
        const amrex::Array4<amrex::Real>& gp = mf_gp.array(mfi);
        const amrex::Array4<const amrex::Real>& gp0 = mf_gp0.const_array(mfi);
        const amrex::Array4<const int>& iblank = mf_iblank.const_array(mfi);
        amrex::ParallelFor(
            gbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                // Replace amr-wind vof values with originals
                if (iblank(i, j, k) == -1) {
                    gp(i, j, k, 0) = gp0(i, j, k, 0);
                    gp(i, j, k, 1) = gp0(i, j, k, 1);
                    gp(i, j, k, 2) = gp0(i, j, k, 2);
                }
            });
    }
}

// Apply pressure gradient to velocity field
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

} // namespace amr_wind::overset_ops

#endif