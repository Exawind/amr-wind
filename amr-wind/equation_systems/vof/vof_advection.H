
#ifndef VOF_ADVECTION_H
#define VOF_ADVECTION_H

#include "amr-wind/equation_systems/vof/vof.H"
#include "amr-wind/equation_systems/vof/SplitAdvection.H"

namespace amr_wind {
namespace pde {

/** Godunov scheme for VOF
 *  \ingroup icns
 */

template <>
struct AdvectionOp<VOF, fvm::Godunov>
{
    AdvectionOp(PDEFields& fields_in, bool has_overset, bool variable_density)
        : fields(fields_in)
        , u_mac(fields_in.repo.get_field("u_mac"))
        , v_mac(fields_in.repo.get_field("v_mac"))
        , w_mac(fields_in.repo.get_field("w_mac"))
    {}

    void operator()(const FieldState fstate, const amrex::Real dt)
    {
        static_assert(
            VOF::ndim == 1, "Invalid number of components for scalar");
        auto& repo = fields.repo;
        auto& geom = repo.mesh().Geom();

        auto& dof_field = fields.field;
        auto bcrec_device = dof_field.bcrec_device();
        //
        // Advect volume using the Explicit Lagrangian onto-cell Sweeping method
        // with PLIC reconstruction
        //

        // Define the sweep time
        isweep += 1;
        if (isweep > 3) isweep = 1;

        for (int lev = 0; lev < repo.num_active_levels(); ++lev) {
            amrex::MFItInfo mfi_info;
            if (amrex::Gpu::notInLaunchRegion())
                mfi_info.EnableTiling(amrex::IntVect(1024, 1024, 1024))
                    .SetDynamic(true);
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
            for (amrex::MFIter mfi(dof_field(lev), mfi_info); mfi.isValid();
                 ++mfi) {
                const auto& bx = mfi.tilebox();
                amrex::FArrayBox tmpfab(amrex::grow(bx, 1), 3 * VOF::ndim);
                tmpfab.setVal(0.0);
                multiphase::split_advection(
                    lev, bx, isweep, dof_field(lev).array(mfi),
                    u_mac(lev).const_array(mfi), v_mac(lev).const_array(mfi),
                    w_mac(lev).const_array(mfi),
                    dof_field.bcrec_device().data(), tmpfab.dataPtr(), geom,
                    dt);

                amrex::Gpu::streamSynchronize();
            }
        }
    }

    PDEFields& fields;
    Field& u_mac;
    Field& v_mac;
    Field& w_mac;
    int isweep = 0;
};

} // namespace pde
} // namespace amr_wind
#endif