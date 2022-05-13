
#ifndef VOF_ADVECTION_H
#define VOF_ADVECTION_H

#include "amr-wind/equation_systems/vof/vof.H"
#include "amr-wind/equation_systems/vof/SplitAdvection.H"
#include "AMReX_MultiFabUtil.H"

namespace amr_wind {
namespace pde {

/** Godunov scheme for VOF
 *  \ingroup vof
 */
template <>
struct AdvectionOp<VOF, fvm::Godunov>
{
    AdvectionOp(
        PDEFields& fields_in, bool /*unused*/, bool /*unused*/, bool /*unused*/)
        : fields(fields_in)
        , u_mac(fields_in.repo.get_field("u_mac"))
        , v_mac(fields_in.repo.get_field("v_mac"))
        , w_mac(fields_in.repo.get_field("w_mac"))
    {
        amrex::ParmParse pp_multiphase("VOF");
        pp_multiphase.query("use_lagrangian", m_use_lagrangian);
        pp_multiphase.query("remove_debris", m_rm_debris);
    }

    void preadvect(const FieldState /*unused*/, const amrex::Real /*unused*/) {}

    void operator()(const FieldState /*unused*/, const amrex::Real dt)
    {
        static_assert(
            VOF::ndim == 1, "Invalid number of components for scalar");

        auto& repo = fields.repo;
        const auto& geom = repo.mesh().Geom();

        auto& aa_x = repo.get_field("advalpha_x");
        auto& aa_y = repo.get_field("advalpha_y");
        auto& aa_z = repo.get_field("advalpha_z");

        auto& dof_field = fields.field;
        //
        // Advect volume using either the Explicit Lagrangian onto-cell or
        // Implicit Eulerian Sweeping method with PLIC reconstruction
        //

        auto flux_x =
            repo.create_scratch_field(1, 0, amr_wind::FieldLoc::XFACE);
        auto flux_y =
            repo.create_scratch_field(1, 0, amr_wind::FieldLoc::YFACE);
        auto flux_z =
            repo.create_scratch_field(1, 0, amr_wind::FieldLoc::ZFACE);

        // Scratch field for fluxC
        auto fluxC = repo.create_scratch_field(1, 0, amr_wind::FieldLoc::CELL);

        // Scratch field for vof field
        auto dof_scratch = repo.create_scratch_field(
            1, fvm::Godunov::nghost_state, amr_wind::FieldLoc::CELL);

        // Define the sweep time
        isweep += 1;
        if (isweep > 3) {
            isweep = 1;
        }

        for (int lev = 0; lev < repo.num_active_levels(); ++lev) {
            // Copy vof field at n to scratch field that will update
            amrex::MultiFab::Copy(
                (*dof_scratch)(lev), dof_field(lev), 0, 0, 1,
                fvm::Godunov::nghost_state);

            amrex::MFItInfo mfi_info;
            if (amrex::Gpu::notInLaunchRegion()) {
                mfi_info.EnableTiling(amrex::IntVect(1024, 1024, 1024))
                    .SetDynamic(true);
            }
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
            for (amrex::MFIter mfi(dof_field(lev), mfi_info); mfi.isValid();
                 ++mfi) {
                const auto& bx = mfi.tilebox();
                amrex::FArrayBox tmpfab(amrex::grow(bx, 1), 2 * VOF::ndim + 2);
                tmpfab.setVal<amrex::RunOn::Device>(0.0);

                // Compression term coefficient
                multiphase::cmask_loop(
                    bx, (*dof_scratch)(lev).array(mfi),
                    (*fluxC)(lev).array(mfi), m_use_lagrangian);

                // First stage of vof advection
                multiphase::split_advection_step(
                    lev, bx, isweep + 0, (*dof_scratch)(lev).array(mfi),
                    (*fluxC)(lev).array(mfi), u_mac(lev).const_array(mfi),
                    v_mac(lev).const_array(mfi), w_mac(lev).const_array(mfi),
                    aa_x(lev).array(mfi), aa_y(lev).array(mfi),
                    aa_z(lev).array(mfi), (*flux_x)(lev).array(mfi),
                    (*flux_y)(lev).array(mfi), (*flux_z)(lev).array(mfi),
                    dof_field.bcrec_device().data(), tmpfab.dataPtr(), geom, dt,
                    m_use_lagrangian);

                amrex::Gpu::streamSynchronize();
            }

            (*dof_scratch)(lev).FillBoundary(geom[lev].periodicity());

            for (amrex::MFIter mfi(dof_field(lev), mfi_info); mfi.isValid();
                 ++mfi) {
                const auto& bx = mfi.tilebox();
                amrex::FArrayBox tmpfab(amrex::grow(bx, 1), 2 * VOF::ndim + 2);
                tmpfab.setVal<amrex::RunOn::Device>(0.0);
                multiphase::split_advection_step(
                    lev, bx, isweep + 1, (*dof_scratch)(lev).array(mfi),
                    (*fluxC)(lev).array(mfi), u_mac(lev).const_array(mfi),
                    v_mac(lev).const_array(mfi), w_mac(lev).const_array(mfi),
                    aa_x(lev).array(mfi), aa_y(lev).array(mfi),
                    aa_z(lev).array(mfi), (*flux_x)(lev).array(mfi),
                    (*flux_y)(lev).array(mfi), (*flux_z)(lev).array(mfi),
                    dof_field.bcrec_device().data(), tmpfab.dataPtr(), geom, dt,
                    m_use_lagrangian);

                amrex::Gpu::streamSynchronize();
            }

            (*dof_scratch)(lev).FillBoundary(geom[lev].periodicity());

            for (amrex::MFIter mfi(dof_field(lev), mfi_info); mfi.isValid();
                 ++mfi) {
                const auto& bx = mfi.tilebox();
                amrex::FArrayBox tmpfab(amrex::grow(bx, 1), 2 * VOF::ndim + 2);
                tmpfab.setVal<amrex::RunOn::Device>(0.0);
                multiphase::split_advection_step(
                    lev, bx, isweep + 2, (*dof_scratch)(lev).array(mfi),
                    (*fluxC)(lev).array(mfi), u_mac(lev).const_array(mfi),
                    v_mac(lev).const_array(mfi), w_mac(lev).const_array(mfi),
                    aa_x(lev).array(mfi), aa_y(lev).array(mfi),
                    aa_z(lev).array(mfi), (*flux_x)(lev).array(mfi),
                    (*flux_y)(lev).array(mfi), (*flux_z)(lev).array(mfi),
                    dof_field.bcrec_device().data(), tmpfab.dataPtr(), geom, dt,
                    m_use_lagrangian);

                amrex::Gpu::streamSynchronize();
            }
        }

        amrex::Vector<amrex::Array<amrex::MultiFab*, AMREX_SPACEDIM>> fluxes(
            repo.num_active_levels());
        for (int lev = 0; lev < repo.num_active_levels(); ++lev) {
            fluxes[lev][0] = &(*flux_x)(lev);
            fluxes[lev][1] = &(*flux_y)(lev);
            fluxes[lev][2] = &(*flux_z)(lev);
        }
        // In order to enforce conservation across coarse-fine boundaries we
        // must be sure to average down the fluxes before we use them
        for (int lev = repo.num_active_levels() - 1; lev > 0; --lev) {
            amrex::IntVect rr =
                geom[lev].Domain().size() / geom[lev - 1].Domain().size();
            amrex::average_down_faces(
                GetArrOfConstPtrs(fluxes[lev]), fluxes[lev - 1], rr,
                geom[lev - 1]);
        }

        for (int lev = 0; lev < repo.num_active_levels(); ++lev) {
            const auto dxinv = geom[lev].InvCellSizeArray();

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
            for (amrex::MFIter mfi(dof_field(lev), amrex::TilingIfNotGPU());
                 mfi.isValid(); ++mfi) {
                const auto& bx = mfi.tilebox();
                auto vof = dof_field(lev).array(mfi);
                auto fx = (*flux_x)(lev).const_array(mfi);
                auto fy = (*flux_y)(lev).const_array(mfi);
                auto fz = (*flux_z)(lev).const_array(mfi);
                // Sum fluxes and contribute directly to dof_field (vof^n)
                amrex::ParallelFor(
                    bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                        vof(i, j, k) +=
                            (fx(i, j, k) - fx(i + 1, j, k)) * dxinv[0] +
                            (fy(i, j, k) - fy(i, j + 1, k)) * dxinv[1] +
                            (fz(i, j, k) - fz(i, j, k + 1)) * dxinv[2];
                        // Do clipping
                        vof(i, j, k) =
                            amrex::max(0.0, amrex::min(1.0, vof(i, j, k)));
                    });

                // Remove debris, if desired
                if (m_rm_debris) {
                    multiphase::debris_loop(bx, dof_field(lev).array(mfi));
                }
            }
        }
    }

    PDEFields& fields;
    Field& u_mac;
    Field& v_mac;
    Field& w_mac;
    int isweep = 0;
    bool m_use_lagrangian{false};
    bool m_rm_debris{true};
};

} // namespace pde
} // namespace amr_wind
#endif
