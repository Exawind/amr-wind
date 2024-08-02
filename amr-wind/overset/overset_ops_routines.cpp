#include "amr-wind/overset/overset_ops_routines.H"

namespace amr_wind::overset_ops {
// Populate approximate signed distance function using vof field
void populate_psi(
    amrex::MultiFab& mf_psi,
    const amrex::MultiFab& mf_vof,
    const amrex::Real i_th,
    const amrex::Real asdf_tiny)
{
    const auto& psi = mf_psi.arrays();
    const auto& vof = mf_vof.const_arrays();
    amrex::ParallelFor(
        mf_psi, mf_psi.n_grow,
        [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
            psi[nbx](i, j, k) = asdf(vof[nbx](i, j, k), i_th, asdf_tiny);
        });
}

// Modify a vof field to not have values that barely differ from 0 or 1
void process_vof(amrex::MultiFab& mf_vof, const amrex::Real vof_tol)
{
    const auto& vof = mf_vof.arrays();
    amrex::ParallelFor(
        mf_vof, mf_vof.n_grow,
        [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
            vof[nbx](i, j, k) =
                vof[nbx](i, j, k) < vof_tol
                    ? 0.0
                    : (vof[nbx](i, j, k) > 1. - vof_tol ? 1.
                                                        : vof[nbx](i, j, k));
        });
}

// Combine overset target vof field with current non-overset vof field
void harmonize_vof(
    amrex::MultiFab& mf_vof_target,
    const amrex::MultiFab& mf_vof_original,
    const amrex::iMultiFab& mf_iblank)
{
    const auto& tg_vof = mf_vof_target.arrays();
    const auto& og_vof = mf_vof_original.const_arrays();
    const auto& iblank = mf_iblank.const_arrays();
    amrex::ParallelFor(
        mf_vof_target,
        [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
            // Replace amr-wind vof values with originals
            if (iblank[nbx](i, j, k) != -1) {
                tg_vof[nbx](i, j, k) = og_vof[nbx](i, j, k);
            }
        });
}

// Populate normal vector with special treatment of overset boundary
void populate_normal_vector(
    amrex::MultiFab& mf_normvec,
    const amrex::MultiFab& mf_vof,
    const amrex::iMultiFab& mf_iblank)
{
    const auto& normvec = mf_normvec.arrays();
    const auto& vof = mf_vof.const_arrays();
    const auto& iblank = mf_iblank.const_arrays();
    // Calculate gradients in each direction with centered diff
    amrex::ParallelFor(
        mf_normvec, mf_normvec.n_grow - amrex::IntVect(1),
        [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
            // Neumann condition across nalu bdy
            int ibdy =
                (iblank[nbx](i, j, k) != iblank[nbx](i - 1, j, k)) ? -1 : 0;
            int jbdy =
                (iblank[nbx](i, j, k) != iblank[nbx](i, j - 1, k)) ? -1 : 0;
            int kbdy =
                (iblank[nbx](i, j, k) != iblank[nbx](i, j, k - 1)) ? -1 : 0;
            // no cell should be isolated such that -1 and 1 are needed
            ibdy =
                (iblank[nbx](i, j, k) != iblank[nbx](i + 1, j, k)) ? +1 : ibdy;
            jbdy =
                (iblank[nbx](i, j, k) != iblank[nbx](i, j + 1, k)) ? +1 : jbdy;
            kbdy =
                (iblank[nbx](i, j, k) != iblank[nbx](i, j, k + 1)) ? +1 : kbdy;
            // Calculate normal
            amrex::Real mx, my, mz, mmag;
            multiphase::youngs_finite_difference_normal_neumann(
                i, j, k, ibdy, jbdy, kbdy, vof[nbx], mx, my, mz);
            // Normalize normal
            mmag = std::sqrt(mx * mx + my * my + mz * mz + 1e-20);
            // Save normal
            normvec[nbx](i, j, k, 0) = mx / mmag;
            normvec[nbx](i, j, k, 1) = my / mmag;
            normvec[nbx](i, j, k, 2) = mz / mmag;
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
    const amrex::MultiFab& mf_gp,
    const amrex::MultiFab& mf_density,
    const amrex::Real Gamma,
    const amrex::Real margin,
    const amrex::Real rho1,
    const amrex::Real rho2)
{
    const auto& fx = mf_fx.arrays();
    const auto& fy = mf_fy.arrays();
    const auto& fz = mf_fz.arrays();
    const auto& vof = mf_vof.const_arrays();
    const auto& tg_vof = mf_target_vof.const_arrays();
    const auto& norm = mf_norm.const_arrays();
    const auto& vel = mf_velocity.const_arrays();
    const auto& gp = mf_gp.const_arrays();
    const auto& rho = mf_density.const_arrays();
    amrex::ParallelFor(
        mf_fx, mf_fx.n_grow,
        [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
            // vof flux
            amrex::Real flux = Gamma * alpha_flux(
                                           i, j, k, 0, margin, vof[nbx],
                                           tg_vof[nbx], norm[nbx]);
            fx[nbx](i, j, k, 0) = flux;
            // density flux
            flux *= (rho1 - rho2);
            fx[nbx](i, j, k, 1) = flux;
            // momentum fluxes (dens flux * face vel)
            amrex::Real uf, vf, wf;
            velocity_face(i, j, k, 0, vof[nbx], vel[nbx], uf, vf, wf);
            fx[nbx](i, j, k, 2) = flux * uf;
            fx[nbx](i, j, k, 3) = flux * vf;
            fx[nbx](i, j, k, 4) = flux * wf;
            // pressure gradient fluxes
            gp_rho_face(i, j, k, 0, vof[nbx], gp[nbx], rho[nbx], uf, vf, wf);
            fx[nbx](i, j, k, 5) = flux * uf;
            fx[nbx](i, j, k, 6) = flux * vf;
            fx[nbx](i, j, k, 7) = flux * wf;
            // Turn "on" all flux faces, later modified in
            // process_fluxes_calc_src
            fx[nbx](i, j, k, 8) = 1.0;
        });
    amrex::ParallelFor(
        mf_fy, mf_fy.n_grow,
        [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
            amrex::Real flux = Gamma * alpha_flux(
                                           i, j, k, 1, margin, vof[nbx],
                                           tg_vof[nbx], norm[nbx]);
            fy[nbx](i, j, k, 0) = flux;
            flux *= (rho1 - rho2);
            fy[nbx](i, j, k, 1) = flux;
            amrex::Real uf, vf, wf;
            velocity_face(i, j, k, 1, vof[nbx], vel[nbx], uf, vf, wf);
            fy[nbx](i, j, k, 2) = flux * uf;
            fy[nbx](i, j, k, 3) = flux * vf;
            fy[nbx](i, j, k, 4) = flux * wf;
            gp_rho_face(i, j, k, 1, vof[nbx], gp[nbx], rho[nbx], uf, vf, wf);
            fy[nbx](i, j, k, 5) = flux * uf;
            fy[nbx](i, j, k, 6) = flux * vf;
            fy[nbx](i, j, k, 7) = flux * wf;
            fy[nbx](i, j, k, 8) = 1.0;
        });
    amrex::ParallelFor(
        mf_fz, mf_fz.n_grow,
        [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
            amrex::Real flux = Gamma * alpha_flux(
                                           i, j, k, 2, margin, vof[nbx],
                                           tg_vof[nbx], norm[nbx]);
            fz[nbx](i, j, k, 0) = flux;
            flux *= (rho1 - rho2);
            fz[nbx](i, j, k, 1) = flux;
            amrex::Real uf, vf, wf;
            velocity_face(i, j, k, 2, vof[nbx], vel[nbx], uf, vf, wf);
            fz[nbx](i, j, k, 2) = flux * uf;
            fz[nbx](i, j, k, 3) = flux * vf;
            fz[nbx](i, j, k, 4) = flux * wf;
            gp_rho_face(i, j, k, 2, vof[nbx], gp[nbx], rho[nbx], uf, vf, wf);
            fz[nbx](i, j, k, 5) = flux * uf;
            fz[nbx](i, j, k, 6) = flux * vf;
            fz[nbx](i, j, k, 7) = flux * wf;
            fz[nbx](i, j, k, 8) = 1.0;
        });
}

// Process reinitialization fluxes - zero non-internal to overset region;
// also calculate pressure source / sink term as a function of fluxes
void process_fluxes_calc_src(
    amrex::MultiFab& mf_fx,
    amrex::MultiFab& mf_fy,
    amrex::MultiFab& mf_fz,
    amrex::MultiFab& mf_psource,
    const amrex::MultiFab& mf_vof,
    const amrex::iMultiFab& mf_iblank)
{
    const auto& fx = mf_fx.arrays();
    const auto& fy = mf_fy.arrays();
    const auto& fz = mf_fz.arrays();
    const auto& sp = mf_psource.arrays();
    const auto& vof = mf_vof.const_arrays();
    const auto& iblank = mf_iblank.const_arrays();
    constexpr amrex::Real tiny = std::numeric_limits<amrex::Real>::epsilon();
    // Zero fluxes based on iblank array
    amrex::ParallelFor(
        mf_fx, mf_fx.n_grow, mf_fx.n_comp,
        [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k, int n) noexcept {
            const bool zero_all =
                (iblank[nbx](i - 1, j, k) + iblank[nbx](i, j, k) > -2);
            fx[nbx](i, j, k, n) *= zero_all ? 0. : 1.;
        });
    amrex::ParallelFor(
        mf_fy, mf_fy.n_grow, mf_fy.n_comp,
        [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k, int n) noexcept {
            const bool zero_all =
                (iblank[nbx](i, j - 1, k) + iblank[nbx](i, j, k) > -2);
            fy[nbx](i, j, k, n) *= zero_all ? 0. : 1.;
        });
    amrex::ParallelFor(
        mf_fz, mf_fz.n_grow, mf_fz.n_comp,
        [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k, int n) noexcept {
            const bool zero_all =
                (iblank[nbx](i, j, k - 1) + iblank[nbx](i, j, k) > -2);
            fz[nbx](i, j, k, n) *= zero_all ? 0. : 1.;
        });
    // With knowledge of fluxes, compute pressure source term
    amrex::ParallelFor(
        mf_psource,
        [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
            sp[nbx](i, j, k) =
                gp_flux_tensor(i, j, k, fx[nbx], fy[nbx], fz[nbx], tiny) &&
                normal_reinit_tensor(
                    i, j, k, fx[nbx], fy[nbx], fz[nbx], vof[nbx], tiny);
        });
}

amrex::Real calculate_pseudo_velocity_scale(
    const amrex::iMultiFab& mf_iblank,
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dx,
    const amrex::Real pvmax)
{
    // The dx of this level should be considered if iblank = -1 here
    // Otherwise, set to max possible value for this mesh (level 0 values)
    return (mf_iblank.min(0) == -1) ? std::min(std::min(dx[0], dx[1]), dx[2])
                                    : pvmax;
}

// Calculate a type of CFL by measuring how much % VOF is being removed per cell
amrex::Real calculate_pseudo_dt_flux(
    const amrex::MultiFab& mf_fx,
    const amrex::MultiFab& mf_fy,
    const amrex::MultiFab& mf_fz,
    const amrex::MultiFab& mf_vof,
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& dx,
    const amrex::Real tol)
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
                    pdt_lim = vof(i, j, k) * dx[2] / fz(i, j, k, 0);
                } else if (fz(i, j, k, 0) < -tol && vof(i, j, k - 1) > tol) {
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

// Apply reinitialization fluxes to modify fields
void apply_fluxes(
    const amrex::MultiFab& mf_fx,
    const amrex::MultiFab& mf_fy,
    const amrex::MultiFab& mf_fz,
    const amrex::MultiFab& mf_psource,
    amrex::MultiFab& mf_vof,
    amrex::MultiFab& mf_dens,
    amrex::MultiFab& mf_vel,
    amrex::MultiFab& mf_gp,
    amrex::MultiFab& mf_pressure,
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dx,
    const amrex::Real ptfac,
    const amrex::Real vof_tol)
{
    const auto& fx = mf_fx.const_arrays();
    const auto& fy = mf_fy.const_arrays();
    const auto& fz = mf_fz.const_arrays();
    const auto& sp = mf_psource.const_arrays();
    const auto& vof = mf_vof.arrays();
    const auto& dens = mf_dens.arrays();
    const auto& vel = mf_vel.arrays();
    const auto& gp = mf_gp.arrays();
    const auto& p = mf_pressure.arrays();

    amrex::ParallelFor(
        mf_vof, [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
            const amrex::Real olddens = dens[nbx](i, j, k);
            vof[nbx](i, j, k) +=
                ptfac *
                ((fx[nbx](i + 1, j, k, 0) - fx[nbx](i, j, k, 0)) / dx[0] +
                 (fy[nbx](i, j + 1, k, 0) - fy[nbx](i, j, k, 0)) / dx[1] +
                 (fz[nbx](i, j, k + 1, 0) - fz[nbx](i, j, k, 0)) / dx[2]);
            dens[nbx](i, j, k) +=
                ptfac *
                ((fx[nbx](i + 1, j, k, 1) - fx[nbx](i, j, k, 1)) / dx[0] +
                 (fy[nbx](i, j + 1, k, 1) - fy[nbx](i, j, k, 1)) / dx[1] +
                 (fz[nbx](i, j, k + 1, 1) - fz[nbx](i, j, k, 1)) / dx[2]);
            vel[nbx](i, j, k, 0) =
                1.0 / dens[nbx](i, j, k) *
                (olddens * vel[nbx](i, j, k, 0) +
                 ptfac *
                     ((fx[nbx](i + 1, j, k, 2) - fx[nbx](i, j, k, 2)) / dx[0] +
                      (fy[nbx](i, j + 1, k, 2) - fy[nbx](i, j, k, 2)) / dx[1] +
                      (fz[nbx](i, j, k + 1, 2) - fz[nbx](i, j, k, 2)) / dx[2]));
            vel[nbx](i, j, k, 1) =
                1.0 / dens[nbx](i, j, k) *
                (olddens * vel[nbx](i, j, k, 1) +
                 ptfac *
                     ((fx[nbx](i + 1, j, k, 3) - fx[nbx](i, j, k, 3)) / dx[0] +
                      (fy[nbx](i, j + 1, k, 3) - fy[nbx](i, j, k, 3)) / dx[1] +
                      (fz[nbx](i, j, k + 1, 3) - fz[nbx](i, j, k, 3)) / dx[2]));
            vel[nbx](i, j, k, 2) =
                1.0 / dens[nbx](i, j, k) *
                (olddens * vel[nbx](i, j, k, 2) +
                 ptfac *
                     ((fx[nbx](i + 1, j, k, 4) - fx[nbx](i, j, k, 4)) / dx[0] +
                      (fy[nbx](i, j + 1, k, 4) - fy[nbx](i, j, k, 4)) / dx[1] +
                      (fz[nbx](i, j, k + 1, 4) - fz[nbx](i, j, k, 4)) / dx[2]));
            gp[nbx](i, j, k, 0) +=
                ptfac *
                ((fx[nbx](i + 1, j, k, 5) - fx[nbx](i, j, k, 5)) / dx[0] +
                 (fy[nbx](i, j + 1, k, 5) - fy[nbx](i, j, k, 5)) / dx[1] +
                 (fz[nbx](i, j, k + 1, 5) - fz[nbx](i, j, k, 5)) / dx[2]);
            gp[nbx](i, j, k, 1) +=
                ptfac *
                ((fx[nbx](i + 1, j, k, 6) - fx[nbx](i, j, k, 6)) / dx[0] +
                 (fy[nbx](i, j + 1, k, 6) - fy[nbx](i, j, k, 6)) / dx[1] +
                 (fz[nbx](i, j, k + 1, 6) - fz[nbx](i, j, k, 6)) / dx[2]);
            gp[nbx](i, j, k, 2) +=
                ptfac *
                ((fx[nbx](i + 1, j, k, 7) - fx[nbx](i, j, k, 7)) / dx[0] +
                 (fy[nbx](i, j + 1, k, 7) - fy[nbx](i, j, k, 7)) / dx[1] +
                 (fz[nbx](i, j, k + 1, 7) - fz[nbx](i, j, k, 7)) / dx[2]);

            // Ensure vof is bounded
            vof[nbx](i, j, k) =
                vof[nbx](i, j, k) < vof_tol
                    ? 0.0
                    : (vof[nbx](i, j, k) > 1. - vof_tol ? 1.
                                                        : vof[nbx](i, j, k));
            // Density bounds are enforced elsewhere
        });
    amrex::ParallelFor(
        mf_pressure,
        [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
            p[nbx](i, j, k) += ptfac * sp[nbx](i, j, k);
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
    const auto& dest = mf_dest.arrays();
    const auto& src = mf_src.const_arrays();
    amrex::ParallelFor(
        mf_dest, [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
            dest[nbx](i, j, k) = std::sqrt(
                src[nbx](i, j, k, 0) * src[nbx](i, j, k, 0) +
                src[nbx](i, j, k, 1) * src[nbx](i, j, k, 1) +
                src[nbx](i, j, k, 2) * src[nbx](i, j, k, 2));
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
    const auto& gp = mf_gp.arrays();
    const auto& rho = mf_density.const_arrays();
    const auto& rho0 = mf_refdens.const_arrays();
    const auto& iblank = mf_iblank.const_arrays();
    amrex::ParallelFor(
        mf_gp, mf_gp.n_grow,
        [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
            if (iblank[nbx](i, j, k) == -1) {
                const amrex::Real dfac =
                    is_pptb ? rho[nbx](i, j, k) - rho0[nbx](i, j, k)
                            : rho[nbx](i, j, k);
                gp[nbx](i, j, k, 0) = 0.;
                gp[nbx](i, j, k, 1) = 0.;
                gp[nbx](i, j, k, 2) = dfac * grav_z;
            }
        });
}

// Swap pressure gradient values in overset region
void replace_gradp(
    amrex::MultiFab& mf_gp,
    const amrex::MultiFab& mf_gp0,
    const amrex::iMultiFab& mf_iblank)
{
    const auto gp = mf_gp.arrays();
    const auto gp0 = mf_gp0.const_arrays();
    const auto iblank = mf_iblank.const_arrays();
    amrex::ParallelFor(
        mf_gp, mf_gp.n_grow,
        [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
            if (iblank[nbx](i, j, k) == -1) {
                gp[nbx](i, j, k, 0) = gp0[nbx](i, j, k, 0);
                gp[nbx](i, j, k, 1) = gp0[nbx](i, j, k, 1);
                gp[nbx](i, j, k, 2) = gp0[nbx](i, j, k, 2);
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
    const auto& vel = mf_vel.arrays();
    const auto& rho = mf_density.const_arrays();
    const auto& gp = mf_gp.const_arrays();
    amrex::ParallelFor(
        mf_vel, [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
            const amrex::Real soverrho = scaling_factor / rho[nbx](i, j, k);
            vel[nbx](i, j, k, 0) -= gp[nbx](i, j, k, 0) * soverrho;
            vel[nbx](i, j, k, 1) -= gp[nbx](i, j, k, 1) * soverrho;
            vel[nbx](i, j, k, 2) -= gp[nbx](i, j, k, 2) * soverrho;
        });
}

} // namespace amr_wind::overset_ops