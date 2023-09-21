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
    if (psi_upw * psi_dnw >= 0.0) {
        // Upwind and downwind are same majority
        if (psi_upw < 0.0) {
            // Upwind in majority gas (avoid over-emptying)
            psi_eval = psi_upw;
        }
        if (psi_dnw > 0.0) {
            // Downwind in majority liquid (avoid over-filling)
            psi_eval = psi_dnw;
        }
    }
    const amrex::Real omtanh2_choose =
        1.0 - tanh(0.5 * psi_eval / i_th) * tanh(0.5 * psi_eval / i_th);
    const amrex::Real a_flux =
        (i_th * gphi - 0.25 * omtanh2_choose * gpsi_norm_avg);
    return a_flux;
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
        const auto& zbx = amrex::surroundingNodes(vbx, 2);
        const amrex::Array4<amrex::Real>& fx = mf_fx.array(mfi);
        const amrex::Array4<amrex::Real>& fy = mf_fy.array(mfi);
        const amrex::Array4<amrex::Real>& fz = mf_fz.array(mfi);
        const amrex::Array4<const amrex::Real>& gpsi = mf_gpsi.const_array(mfi);
        const amrex::Array4<const amrex::Real>& vof = mf_vof.const_array(mfi);
        // Populate vof and density fluxes for each direction
        amrex::ParallelFor(
            xbx, ybx, zbx,
            [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                // vof flux
                amrex::Real f_vof = alpha_flux(i, j, k, 0, i_th, vof, gpsi, dx);
                fx(i, j, k, 0) = f_vof;
                // density flux
                // fx(i, j, k, 1) = rho1 * fx(i,j,k,0) + rho2 * (1.0 -
                // fx(i,j,k,0)); pressure gradient fluxes
            },
            [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                amrex::Real f_vof = alpha_flux(i, j, k, 1, i_th, vof, gpsi, dx);
                fy(i, j, k, 0) = f_vof;
            },
            [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                amrex::Real f_vof = alpha_flux(i, j, k, 2, i_th, vof, gpsi, dx);
                fz(i, j, k, 0) = f_vof;
            });
    }
}

void process_fluxes(
    amrex::MultiFab& mf_fx,
    amrex::MultiFab& mf_fy,
    amrex::MultiFab& mf_fz,
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
        const amrex::Array4<const int>& iblank = mf_iblank.const_array(mfi);
        // Only faces with iblank = -1 on both sides can have nonzero flux
        // Faces with single-phase on one side cannot have flux
        amrex::ParallelFor(
            xbx, ybx, zbx,
            [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                bool zero_all = (iblank(i - 1, j, k) + iblank(i, j, k) > -2);
                for (int n = 0; n < ncomp; ++n) {
                    fx(i, j, k, n) *= zero_all ? 0. : 1.;
                }
            },
            [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                bool zero_all = (iblank(i, j - 1, k) + iblank(i, j, k) > -2);
                for (int n = 0; n < ncomp; ++n) {
                    fy(i, j, k, n) *= zero_all ? 0. : 1.;
                }
            },
            [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                bool zero_all = (iblank(i, j, k - 1) + iblank(i, j, k) > -2);
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
        // Only faces with iblank = -1 on both sides can have nonzero flux
        amrex::ParallelFor(
            vbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                const amrex::Real oldvof = vof(i, j, k);
                vof(i, j, k) +=
                    i_th * ((fx(i + 1, j, k, 0) - fx(i, j, k, 0)) / dx[0] +
                            (fy(i, j + 1, k, 0) - fy(i, j, k, 0)) / dx[1] +
                            (fz(i, j, k + 1, 0) - fz(i, j, k, 0)) / dx[2]);
                /*if (i == 20 && j == 20 && vof(i, j, k) > 1.0) {
                    std::cout << i << " " << j << " " << k << " "
                              << vof(i, j, k) << " " << oldvof << " "
                              << fz(i, j, k + 1, 0) << " " << fz(i, j, k, 0)
                              << std::endl;
                }*/
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
            const amrex::Real i_th = 0.125 * std::cbrt(dx[0] * dx[1] * dx[2]);

            // Populate approximate signed distance function
            populate_grad_psi(
                (*grad_asdf)(lev), vof(lev), iblank_cell(lev), dx, i_th);

            // Sharpening fluxes for vof, density, pressure gradient
            populate_sharpen_fluxes(
                (*flux_x)(lev), (*flux_y)(lev), (*flux_z)(lev),
                (*grad_asdf)(lev), vof(lev), dx, i_th, rhol, rhog);

            // Zero fluxes outside nalu portion of domain
            process_fluxes(
                (*flux_x)(lev), (*flux_y)(lev), (*flux_z)(lev),
                iblank_cell(lev));

            // Apply fluxes
            apply_fluxes(
                (*flux_x)(lev), (*flux_y)(lev), (*flux_z)(lev), vof(lev), dx,
                i_th);

            // Replace vof with another field to look at it
            // equate_vof((*grad_asdf)(lev), vof(lev));
        }

        // Average down fluxes
        // Get pseudo CFL (multiplier for sharpening term)
        // Apply fluxes

        // Average down vof for consistency

        // Fillpatch for ghost cells
        vof.fillpatch(sim.time().current_time());

        // Update density
    }
}
} // namespace overset
} // namespace amr_wind

#endif