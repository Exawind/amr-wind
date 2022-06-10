
#ifndef VOF_ADVECTION_H
#define VOF_ADVECTION_H

#include "amr-wind/equation_systems/vof/vof.H"
#include "amr-wind/equation_systems/vof/SplitAdvection.H"

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

        // cppcheck-suppress constVariable
        auto& dof_field = fields.field;
        //
        // Advect volume using either the Explicit Lagrangian onto-cell or
        // Implicit Eulerian Sweeping method with PLIC reconstruction
        //

        // Scratch field for fluxC
        auto fluxC = repo.create_scratch_field(1, 0, amr_wind::FieldLoc::CELL);

        // Define the sweep time
        isweep += 1;
        if (isweep > 3) {
            isweep = 1;
        }

        for (int lev = 0; lev < repo.num_active_levels(); ++lev) {
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
                amrex::FArrayBox tmpfab(amrex::grow(bx, 1), 2 * VOF::ndim);
                tmpfab.setVal<amrex::RunOn::Device>(0.0);

                multiphase::cmask_loop(
                    bx, dof_field(lev).array(mfi), (*fluxC)(lev).array(mfi),
                    m_use_lagrangian);

                multiphase::split_advection_step(
                    lev, bx, isweep + 0, dof_field(lev).array(mfi),
                    (*fluxC)(lev).array(mfi), u_mac(lev).const_array(mfi),
                    v_mac(lev).const_array(mfi), w_mac(lev).const_array(mfi),
                    dof_field.bcrec_device().data(), tmpfab.dataPtr(), geom, dt,
                    m_use_lagrangian);

                amrex::Gpu::streamSynchronize();
            }

            dof_field(lev).FillBoundary(geom[lev].periodicity());

            for (amrex::MFIter mfi(dof_field(lev), mfi_info); mfi.isValid();
                 ++mfi) {
                const auto& bx = mfi.tilebox();
                amrex::FArrayBox tmpfab(amrex::grow(bx, 1), 2 * VOF::ndim);
                tmpfab.setVal<amrex::RunOn::Device>(0.0);
                multiphase::split_advection_step(
                    lev, bx, isweep + 1, dof_field(lev).array(mfi),
                    (*fluxC)(lev).array(mfi), u_mac(lev).const_array(mfi),
                    v_mac(lev).const_array(mfi), w_mac(lev).const_array(mfi),
                    dof_field.bcrec_device().data(), tmpfab.dataPtr(), geom, dt,
                    m_use_lagrangian);

                amrex::Gpu::streamSynchronize();
            }

            dof_field(lev).FillBoundary(geom[lev].periodicity());

            for (amrex::MFIter mfi(dof_field(lev), mfi_info); mfi.isValid();
                 ++mfi) {
                const auto& bx = mfi.tilebox();
                amrex::FArrayBox tmpfab(amrex::grow(bx, 1), 2 * VOF::ndim);
                tmpfab.setVal<amrex::RunOn::Device>(0.0);
                multiphase::split_advection_step(
                    lev, bx, isweep + 2, dof_field(lev).array(mfi),
                    (*fluxC)(lev).array(mfi), u_mac(lev).const_array(mfi),
                    v_mac(lev).const_array(mfi), w_mac(lev).const_array(mfi),
                    dof_field.bcrec_device().data(), tmpfab.dataPtr(), geom, dt,
                    m_use_lagrangian);

                if (m_rm_debris) {
                    multiphase::debris_loop(bx, dof_field(lev).array(mfi));
                }
                amrex::Gpu::streamSynchronize();
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
