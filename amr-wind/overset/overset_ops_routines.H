#ifndef OVERSET_OPS_ROUTINES_H_
#define OVERSET_OPS_ROUTINES_H_

#include "amr-wind/equation_systems/vof/volume_fractions.H"
#include "amr-wind/overset/overset_ops_K.H"

namespace amr_wind::overset_ops {

// Populate approximate signed distance function using vof field
void populate_psi(
    amrex::MultiFab& mf_psi,
    const amrex::MultiFab& mf_vof,
    const amrex::Real i_th,
    const amrex::Real asdf_tiny)
{
    const auto& psi_arrs = mf_psi.arrays();
    const auto& vof_arrs = mf_vof.const_arrays();
    amrex::ParallelFor(
        mf_psi, mf_psi.n_grow,
        [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
            psi_arrs[nbx](i, j, k) =
                asdf(vof_arrs[nbx](i, j, k), i_th, asdf_tiny);
        });
}

// Modify a vof field to not have values that barely differ from 0 or 1
void process_vof(amrex::MultiFab& mf_vof, const amrex::Real vof_tol)
{
    const auto& vof_arrs = mf_vof.arrays();
    amrex::ParallelFor(
        mf_vof, mf_vof.n_grow,
        [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
            vof_arrs[nbx](i, j, k) =
                vof_arrs[nbx](i, j, k) < vof_tol
                    ? 0.0
                    : (vof_arrs[nbx](i, j, k) > 1. - vof_tol
                           ? 1.
                           : vof_arrs[nbx](i, j, k));
        });
}

// Combine overset target vof field with current non-overset vof field
void harmonize_vof(
    amrex::MultiFab& mf_vof_target,
    const amrex::MultiFab& mf_vof_original,
    const amrex::iMultiFab& mf_iblank)
{
    const auto& tg_vof_arrs = mf_vof_target.arrays();
    const auto& og_vof_arrs = mf_vof_original.const_arrays();
    const auto& iblank_arrs = mf_iblank.const_arrays();
    amrex::ParallelFor(
        mf_vof_target,
        [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
            // Replace amr-wind vof values with originals
            if (iblank_arrs[nbx](i, j, k) != -1) {
                tg_vof_arrs[nbx](i, j, k) = og_vof_arrs[nbx](i, j, k);
            }
        });
}

// Populate normal vector with special treatment of overset boundary
void populate_normal_vector(
    amrex::MultiFab& mf_normvec,
    const amrex::MultiFab& mf_vof,
    const amrex::iMultiFab& mf_iblank)
{
    const auto& normvec_arrs = mf_normvec.arrays();
    const auto& vof_arrs = mf_vof.const_arrays();
    const auto& iblank_arrs = mf_iblank.const_arrays();
    // Calculate gradients in each direction with centered diff
    amrex::ParallelFor(
        mf_normvec, mf_normvec.n_grow - amrex::IntVect(1),
        [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
            // Neumann condition across nalu bdy
            int ibdy =
                (iblank_arrs[nbx](i, j, k) != iblank_arrs[nbx](i - 1, j, k))
                    ? -1
                    : 0;
            int jbdy =
                (iblank_arrs[nbx](i, j, k) != iblank_arrs[nbx](i, j - 1, k))
                    ? -1
                    : 0;
            int kbdy =
                (iblank_arrs[nbx](i, j, k) != iblank_arrs[nbx](i, j, k - 1))
                    ? -1
                    : 0;
            // no cell should be isolated such that -1 and 1 are needed
            ibdy = (iblank_arrs[nbx](i, j, k) != iblank_arrs[nbx](i + 1, j, k))
                       ? +1
                       : ibdy;
            jbdy = (iblank_arrs[nbx](i, j, k) != iblank_arrs[nbx](i, j + 1, k))
                       ? +1
                       : jbdy;
            kbdy = (iblank_arrs[nbx](i, j, k) != iblank_arrs[nbx](i, j, k + 1))
                       ? +1
                       : kbdy;
            // Calculate normal
            amrex::Real mx, my, mz, mmag;
            multiphase::youngs_finite_difference_normal_neumann(
                i, j, k, ibdy, jbdy, kbdy, vof_arrs[nbx], mx, my, mz);
            // Normalize normal
            mmag = std::sqrt(mx * mx + my * my + mz * mz + 1e-20);
            // Save normal
            normvec_arrs[nbx](i, j, k, 0) = mx / mmag;
            normvec_arrs[nbx](i, j, k, 1) = my / mmag;
            normvec_arrs[nbx](i, j, k, 2) = mz / mmag;
        });
}

// Calculate fluxes for reinitialization over entire domain without concern for
// overset bdy
void populate_sharpen_fluxes(
    amrex::MultiFab& mf_fx,
    amrex::MultiFab& mf_fy,
    amrex::MultiFab& mf_fz,
    const amrex::MultiFab& mf_vof,
    const amrex::MultiFab& mf_target_vof,
    const amrex::MultiFab& mf_norm,
    const amrex::MultiFab& mf_velocity,
    const amrex::Real margin,
    const amrex::Real rho1,
    const amrex::Real rho2)
{
    for (amrex::MFIter mfi(mf_vof); mfi.isValid(); ++mfi) {
        const auto& vbx = mfi.validbox();
        const auto& xbx = amrex::surroundingNodes(vbx, 0);
        const auto& ybx = amrex::surroundingNodes(vbx, 1);
        const auto& zbx = amrex::surroundingNodes(vbx, 2);
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
            xbx, ybx, zbx,
            [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                // vof flux
                amrex::Real flux =
                    alpha_flux(i, j, k, 0, margin, vof, tg_vof, norm);
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
            },
            [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                amrex::Real flux =
                    alpha_flux(i, j, k, 1, margin, vof, tg_vof, norm);
                fy(i, j, k, 0) = flux;
                flux *= (rho1 - rho2);
                fy(i, j, k, 1) = flux;
                amrex::Real uf, vf, wf;
                velocity_face(i, j, k, 1, vof, vel, uf, vf, wf);
                fy(i, j, k, 2) = flux * uf;
                fy(i, j, k, 3) = flux * vf;
                fy(i, j, k, 4) = flux * wf;
            },
            [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                amrex::Real flux =
                    alpha_flux(i, j, k, 2, margin, vof, tg_vof, norm);
                fz(i, j, k, 0) = flux;
                flux *= (rho1 - rho2);
                fz(i, j, k, 1) = flux;
                amrex::Real uf, vf, wf;
                velocity_face(i, j, k, 2, vof, vel, uf, vf, wf);
                fz(i, j, k, 2) = flux * uf;
                fz(i, j, k, 3) = flux * vf;
                fz(i, j, k, 4) = flux * wf;
            });
    }
}

// Process reinitialization fluxes - zero non-internal to overset region
void process_fluxes(
    amrex::MultiFab& mf_fx,
    amrex::MultiFab& mf_fy,
    amrex::MultiFab& mf_fz,
    const amrex::iMultiFab& mf_iblank)
{
    int ncompx = mf_fx.nComp();
    int ncompy = mf_fy.nComp();
    int ncompz = mf_fz.nComp();
    for (amrex::MFIter mfi(mf_iblank); mfi.isValid(); ++mfi) {
        const auto& vbx = mfi.validbox();
        const auto& xbx = amrex::surroundingNodes(vbx, 0);
        const auto& ybx = amrex::surroundingNodes(vbx, 1);
        const auto& zbx = amrex::surroundingNodes(vbx, 2);
        const amrex::Array4<amrex::Real>& fx = mf_fx.array(mfi);
        const amrex::Array4<amrex::Real>& fy = mf_fy.array(mfi);
        const amrex::Array4<amrex::Real>& fz = mf_fz.array(mfi);
        const amrex::Array4<const int>& iblank = mf_iblank.const_array(mfi);
        // Only faces with iblank = -1 on both sides can have nonzero flux
        amrex::ParallelFor(
            xbx, ybx, zbx,
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
    }
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
    amrex::MultiFab& mf_vof,
    amrex::MultiFab& mf_dens,
    amrex::MultiFab& mf_vel,
    amrex::Real ptfac,
    amrex::Real vof_tol)
{
    const auto& fx = mf_fx.const_arrays();
    const auto& fy = mf_fy.const_arrays();
    const auto& fz = mf_fz.const_arrays();
    const auto& vof = mf_vof.arrays();
    const auto& dens = mf_dens.arrays();
    const auto& vel = mf_vel.arrays();

    amrex::ParallelFor(
        mf_vof, [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
            const amrex::Real olddens = dens[nbx](i, j, k);
            vof[nbx](i, j, k) +=
                ptfac * ((fx[nbx](i + 1, j, k, 0) - fx[nbx](i, j, k, 0)) +
                         (fy[nbx](i, j + 1, k, 0) - fy[nbx](i, j, k, 0)) +
                         (fz[nbx](i, j, k + 1, 0) - fz[nbx](i, j, k, 0)));
            dens[nbx](i, j, k) +=
                ptfac * ((fx[nbx](i + 1, j, k, 1) - fx[nbx](i, j, k, 1)) +
                         (fy[nbx](i, j + 1, k, 1) - fy[nbx](i, j, k, 1)) +
                         (fz[nbx](i, j, k + 1, 1) - fz[nbx](i, j, k, 1)));
            vel[nbx](i, j, k, 0) =
                1.0 / dens[nbx](i, j, k) *
                (olddens * vel[nbx](i, j, k, 0) +
                 ptfac * ((fx[nbx](i + 1, j, k, 2) - fx[nbx](i, j, k, 2)) +
                          (fy[nbx](i, j + 1, k, 2) - fy[nbx](i, j, k, 2)) +
                          (fz[nbx](i, j, k + 1, 2) - fz[nbx](i, j, k, 2))));
            vel[nbx](i, j, k, 1) =
                1.0 / dens[nbx](i, j, k) *
                (olddens * vel[nbx](i, j, k, 1) +
                 ptfac * ((fx[nbx](i + 1, j, k, 3) - fx[nbx](i, j, k, 3)) +
                          (fy[nbx](i, j + 1, k, 3) - fy[nbx](i, j, k, 3)) +
                          (fz[nbx](i, j, k + 1, 3) - fz[nbx](i, j, k, 3))));
            vel[nbx](i, j, k, 2) =
                1.0 / dens[nbx](i, j, k) *
                (olddens * vel[nbx](i, j, k, 2) +
                 ptfac * ((fx[nbx](i + 1, j, k, 4) - fx[nbx](i, j, k, 4)) +
                          (fy[nbx](i, j + 1, k, 4) - fy[nbx](i, j, k, 4)) +
                          (fz[nbx](i, j, k + 1, 4) - fz[nbx](i, j, k, 4))));

            // Ensure vof is bounded
            vof[nbx](i, j, k) =
                vof[nbx](i, j, k) < vof_tol
                    ? 0.0
                    : (vof[nbx](i, j, k) > 1. - vof_tol ? 1.
                                                        : vof[nbx](i, j, k));
            // Density bounds are enforced elsewhere
        });
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
void equate_field(amrex::MultiFab& mf_dest, const amrex::MultiFab& mf_src)
{
    const auto& dest_arrs = mf_dest.arrays();
    const auto& src_arrs = mf_src.const_arrays();
    amrex::ParallelFor(
        mf_dest, [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
            dest_arrs[nbx](i, j, k) = std::sqrt(
                src_arrs[nbx](i, j, k, 0) * src_arrs[nbx](i, j, k, 0) +
                src_arrs[nbx](i, j, k, 1) * src_arrs[nbx](i, j, k, 1) +
                src_arrs[nbx](i, j, k, 2) * src_arrs[nbx](i, j, k, 2));
        });
}

// Replace pressure gradient with hydrostatic field in overset regions
void replace_gradp_hydrostatic(
    amrex::MultiFab& mf_gp,
    const amrex::MultiFab& mf_density,
    const amrex::MultiFab& mf_refdens,
    const amrex::iMultiFab& mf_iblank,
    const amrex::Real grav_z,
    const bool is_pptb)
{
    const auto& gp_arrs = mf_gp.arrays();
    const auto& rho_arrs = mf_density.const_arrays();
    const auto& rho0_arrs = mf_refdens.const_arrays();
    const auto& iblank_arrs = mf_iblank.const_arrays();
    amrex::ParallelFor(
        mf_gp, mf_gp.n_grow,
        [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
            if (iblank_arrs[nbx](i, j, k) == -1) {
                const amrex::Real dfac =
                    is_pptb ? rho_arrs[nbx](i, j, k) - rho0_arrs[nbx](i, j, k)
                            : rho_arrs[nbx](i, j, k);
                gp_arrs[nbx](i, j, k, 0) = 0.;
                gp_arrs[nbx](i, j, k, 1) = 0.;
                gp_arrs[nbx](i, j, k, 2) = dfac * grav_z;
            }
        });
}

// Swap pressure gradient values in overset region
void replace_gradp(
    amrex::MultiFab& mf_gp,
    const amrex::MultiFab& mf_gp0,
    const amrex::iMultiFab& mf_iblank)
{
    const auto gp_arrs = mf_gp.arrays();
    const auto gp0_arrs = mf_gp0.const_arrays();
    const auto iblank_arrs = mf_iblank.const_arrays();
    amrex::ParallelFor(
        mf_gp, mf_gp.n_grow,
        [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
            if (iblank_arrs[nbx](i, j, k) == -1) {
                gp_arrs[nbx](i, j, k, 0) = gp0_arrs[nbx](i, j, k, 0);
                gp_arrs[nbx](i, j, k, 1) = gp0_arrs[nbx](i, j, k, 1);
                gp_arrs[nbx](i, j, k, 2) = gp0_arrs[nbx](i, j, k, 2);
            }
        });
}

// Apply pressure gradient to velocity field
void apply_pressure_gradient(
    amrex::MultiFab& mf_vel,
    const amrex::MultiFab& mf_density,
    const amrex::MultiFab& mf_gp,
    const amrex::Real scaling_factor)
{
    const auto& vel_arrs = mf_vel.arrays();
    const auto& rho_arrs = mf_density.const_arrays();
    const auto& gp_arrs = mf_gp.const_arrays();
    amrex::ParallelFor(
        mf_vel, [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
            const amrex::Real soverrho =
                scaling_factor / rho_arrs[nbx](i, j, k);
            vel_arrs[nbx](i, j, k, 0) -= gp_arrs[nbx](i, j, k, 0) * soverrho;
            vel_arrs[nbx](i, j, k, 1) -= gp_arrs[nbx](i, j, k, 1) * soverrho;
            vel_arrs[nbx](i, j, k, 2) -= gp_arrs[nbx](i, j, k, 2) * soverrho;
        });
}

} // namespace amr_wind::overset_ops

#endif