#ifndef VOF_MOMENTUM_FLUX_H
#define VOF_MOMENTUM_FLUX_H

#include "amr-wind/equation_systems/vof/volume_fractions.H"

namespace amr_wind::multiphase {

static void hybrid_fluxes(
    const FieldRepo& repo,
    const int ncomp,
    const amrex::Gpu::DeviceVector<int> iconserv,
    ScratchField& flux_x,
    ScratchField& flux_y,
    ScratchField& flux_z,
    const Field& dof_field,
    const Field& src_term,
    const Field& rho_o,
    const Field& rho_nph,
    const Field& u_mac,
    const Field& v_mac,
    const Field& w_mac,
    amrex::BCRec const* velbc,
    amrex::BCRec const* rhobc,
    const amrex::Real dt,
    godunov::scheme mflux_scheme)
{
    // Get geometry
    const auto& geom = repo.mesh().Geom();
    // Get advected alpha fields
    const auto& advalpha_x = repo.get_field("advalpha_x");
    const auto& advalpha_y = repo.get_field("advalpha_y");
    const auto& advalpha_z = repo.get_field("advalpha_z");
    // Get VOF
    const auto& vof = repo.get_field("vof");
    // Get phase densities (how?)
    amrex::Real rho1 = 1000.0;
    amrex::Real rho2 = 1.0;

    // Create scratch arrays for local flux storage
    auto ftmp_x =
        repo.create_scratch_field(ncomp, 0, amr_wind::FieldLoc::XFACE);
    auto ftmp_y =
        repo.create_scratch_field(ncomp, 0, amr_wind::FieldLoc::YFACE);
    auto ftmp_z =
        repo.create_scratch_field(ncomp, 0, amr_wind::FieldLoc::ZFACE);

    // Create iconserv = 0 array for density, avoid multiplication by face vel
    amrex::Gpu::DeviceVector<int> idnsty;
    idnsty.resize(1, 0);

    for (int lev = 0; lev < repo.num_active_levels(); ++lev) {
        // Set up temporary arrays for fluxes that will replace other quantities
        // form multifab for transport variable and source term
        amrex::MultiFab q(
            dof_field(lev).boxArray(), dof_field(lev).DistributionMap(), ncomp,
            fvm::Godunov::nghost_state);
        amrex::MultiFab::Copy(
            q, dof_field(lev), 0, 0, ncomp, fvm::Godunov::nghost_state);
        amrex::MultiFab fq(
            src_term(lev).boxArray(), src_term(lev).DistributionMap(), ncomp,
            fvm::Godunov::nghost_src);
        amrex::MultiFab::Copy(
            fq, src_term(lev), 0, 0, ncomp, fvm::Godunov::nghost_src);
        amrex::MultiFab frho(
            src_term(lev).boxArray(), src_term(lev).DistributionMap(), 1,
            fvm::Godunov::nghost_src);

        amrex::MultiFab::Multiply(
            q, rho_o(lev), 0, 0, 1, fvm::Godunov::nghost_state);
        amrex::MultiFab::Multiply(
            q, rho_o(lev), 0, 1, 1, fvm::Godunov::nghost_state);
        amrex::MultiFab::Multiply(
            q, rho_o(lev), 0, 2, 1, fvm::Godunov::nghost_state);
        amrex::MultiFab::Multiply(
            fq, rho_nph(lev), 0, 0, 1, fvm::Godunov::nghost_src);
        amrex::MultiFab::Multiply(
            fq, rho_nph(lev), 0, 1, 1, fvm::Godunov::nghost_src);
        amrex::MultiFab::Multiply(
            fq, rho_nph(lev), 0, 2, 1, fvm::Godunov::nghost_src);
        frho.setVal(0.0);

        amrex::MFItInfo mfi_info;
        if (amrex::Gpu::notInLaunchRegion()) {
            mfi_info.EnableTiling(amrex::IntVect(1024, 1024, 1024))
                .SetDynamic(true);
        }
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        for (amrex::MFIter mfi(dof_field(lev), mfi_info); mfi.isValid();
             ++mfi) {
            const auto& bx = mfi.tilebox();
            const auto& bxg1 = amrex::grow(bx, 1);
            const auto& xbx = amrex::surroundingNodes(bx, 0);
            const auto& ybx = amrex::surroundingNodes(bx, 1);
            const auto& zbx = amrex::surroundingNodes(bx, 2);
            // Memory for working arrays in godunov::compute_fluxes
            amrex::FArrayBox tmpfab(amrex::grow(bx, 1), ncomp * 14);
            // Arrays for easy reference (w for working array)
            auto F_x = flux_x(lev).array(mfi);
            auto F_y = flux_y(lev).array(mfi);
            auto F_z = flux_z(lev).array(mfi);
            auto Fw_x = (*ftmp_x)(lev).array(mfi);
            auto Fw_y = (*ftmp_y)(lev).array(mfi);
            auto Fw_z = (*ftmp_z)(lev).array(mfi);
            auto aa_x = advalpha_x(lev).const_array(mfi);
            auto aa_y = advalpha_y(lev).const_array(mfi);
            auto aa_z = advalpha_z(lev).const_array(mfi);
            // Compute momentum flux quantities (multiplied by face velocity)
            godunov::compute_fluxes(
                lev, bx, ncomp, Fw_x, Fw_y, Fw_z, q.const_array(mfi),
                u_mac(lev).const_array(mfi), v_mac(lev).const_array(mfi),
                w_mac(lev).const_array(mfi), fq.const_array(mfi), velbc,
                iconserv.data(), tmpfab.dataPtr(), geom, dt, mflux_scheme);

            // Where interface is present, replace current momentum flux
            // quantities with those from mflux_scheme
            auto volfrac = vof(lev).const_array(mfi);
            amrex::ParallelFor(
                bxg1, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    bool vf_i = multiphase::interface_band(i, j, k, volfrac);
                    // Fluxes in each direction, check with flag too
                    if (xbx.contains(i, j, k)) {
                        bool vf_nb =
                            multiphase::interface_band(i - 1, j, k, volfrac);
                        if (vf_i || vf_nb) {
                            for (int n = 0; n < ncomp; ++n) {
                                F_x(i, j, k, n) = Fw_x(i, j, k, n);
                            }
                        }
                    }
                    if (ybx.contains(i, j, k)) {
                        bool vf_nb =
                            multiphase::interface_band(i, j - 1, k, volfrac);
                        if (vf_i || vf_nb) {
                            for (int n = 0; n < ncomp; ++n) {
                                F_y(i, j, k, n) = Fw_y(i, j, k, n);
                            }
                        }
                    }
                    if (zbx.contains(i, j, k)) {
                        bool vf_nb =
                            multiphase::interface_band(i, j, k - 1, volfrac);
                        if (vf_i || vf_nb) {
                            for (int n = 0; n < ncomp; ++n) {
                                F_z(i, j, k, n) = Fw_z(i, j, k, n);
                            }
                        }
                    }
                });

            // Compute density flux quantities (not multiplied by face velocity)
            godunov::compute_fluxes(
                lev, bx, 1, Fw_x, Fw_y, Fw_z, rho_o(lev).const_array(mfi),
                u_mac(lev).const_array(mfi), v_mac(lev).const_array(mfi),
                w_mac(lev).const_array(mfi), frho.const_array(mfi), rhobc,
                idnsty.data(), tmpfab.dataPtr(), geom, dt, mflux_scheme);

            // When interface is present, multiply by consistent density and
            // divide by interpolated density to replace flux value with
            // mass-consistent, Favre-averaged flux
            amrex::ParallelFor(
                bxg1, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    bool vf_i = multiphase::interface_band(i, j, k, volfrac);
                    // Fluxes in each direction, check with flag too
                    if (xbx.contains(i, j, k)) {
                        bool vf_nb =
                            multiphase::interface_band(i - 1, j, k, volfrac);
                        if (vf_i || vf_nb) {
                            amrex::Real rhof = rho1 * aa_x(i, j, k) +
                                               rho2 * (1.0 - aa_x(i, j, k));
                            for (int n = 0; n < ncomp; ++n) {
                                F_x(i, j, k, n) *= rhof / Fw_x(i, j, k, 0);
                            }
                        }
                    }
                    if (ybx.contains(i, j, k)) {
                        bool vf_nb =
                            multiphase::interface_band(i, j - 1, k, volfrac);
                        if (vf_i || vf_nb) {
                            amrex::Real rhof = rho1 * aa_y(i, j, k) +
                                               rho2 * (1.0 - aa_y(i, j, k));
                            for (int n = 0; n < ncomp; ++n) {
                                F_y(i, j, k, n) *= rhof / Fw_y(i, j, k, 0);
                            }
                        }
                    }
                    if (zbx.contains(i, j, k)) {
                        bool vf_nb =
                            multiphase::interface_band(i, j, k - 1, volfrac);
                        if (vf_i || vf_nb) {
                            amrex::Real rhof = rho1 * aa_z(i, j, k) +
                                               rho2 * (1.0 - aa_z(i, j, k));
                            for (int n = 0; n < ncomp; ++n) {
                                F_z(i, j, k, n) *= rhof / Fw_z(i, j, k, 0);
                            }
                        }
                    }
                });
        }
    }
}
} // namespace amr_wind::multiphase

#endif
