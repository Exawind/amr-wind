
#ifndef VOF_ADVECTION_H
#define VOF_ADVECTION_H

#include "amr-wind/equation_systems/vof/vof.H"
#include "amr-wind/equation_systems/vof/SplitLagrangianAdvection.H"

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
        auto& repo = fields.repo;
        auto& geom = repo.mesh().Geom();
        auto& dof_field = fields.field.state(fstate);
        auto bcrec_device = dof_field.bcrec_device();
        //
        // Advect volume using the Explicit Lagrangian onto-cell Sweeping method
        // with PLIC reconstruction
        //
        for (int lev = 0; lev < repo.num_active_levels(); ++lev) {
            u_mac(lev).FillBoundary(geom[lev].periodicity());
            v_mac(lev).FillBoundary(geom[lev].periodicity());
            w_mac(lev).FillBoundary(geom[lev].periodicity());

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

                multiphase::split_lagrangian_advection(
                    lev, bx, dof_field(lev).array(mfi),
                    u_mac(lev).const_array(mfi), v_mac(lev).const_array(mfi),
                    w_mac(lev).const_array(mfi),
                    dof_field.bcrec_device().data(), geom, dt);

                amrex::Gpu::streamSynchronize();
            }
        }
    }

    PDEFields& fields;
    Field& u_mac;
    Field& v_mac;
    Field& w_mac;
};

} // namespace pde
} // namespace amr_wind
#endif