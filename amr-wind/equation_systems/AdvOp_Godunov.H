#ifndef ADVOP_GODUNOV_H
#define ADVOP_GODUNOV_H

#include <type_traits>

#include "amr-wind/convection/Godunov.H"
#include "amr-wind/equation_systems/SchemeTraits.H"
#include "amr-wind/equation_systems/PDETraits.H"
#include "amr-wind/equation_systems/PDEOps.H"

#include "AMReX_Gpu.H"
#include "AMReX_ParmParse.H"

namespace amr_wind {
namespace pde {

/** Godunov advection operator for scalar transport equations
 *  \ingroup pdeop
 */
template <typename PDE>
struct AdvectionOp<
    PDE,
    fvm::Godunov,
    typename std::enable_if<std::is_base_of<ScalarTransport, PDE>::value>::type>
{
    AdvectionOp(PDEFields& fields_in)
    : fields(fields_in)
    , density(fields_in.repo.get_field("density"))
    , u_mac(fields_in.repo.get_field("u_mac"))
    , v_mac(fields_in.repo.get_field("v_mac"))
    , w_mac(fields_in.repo.get_field("w_mac"))
    {
        amrex::ParmParse pp("incflo");
        pp.query("godunov_ppm", godunov_ppm);

        // TODO: Need iconserv flag to be adjusted???
        iconserv.resize(PDE::ndim, 1);
    }

    void operator()(const FieldState fstate, const amrex::Real dt,
                    const bool /* has_overset */)
    {
        static_assert(PDE::ndim == 1, "Invalid number of components for scalar");
        auto& repo = fields.repo;
        auto& geom = repo.mesh().Geom();

        auto& src_term = fields.src_term;
        auto& conv_term = fields.conv_term;
        auto& dof_field = fields.field.state(fstate);

        // only needed if multiplying by rho below
        auto& den = density.state(fstate);

        for (int lev = 0; lev < repo.num_active_levels(); ++lev) {
            amrex::MFItInfo mfi_info;
            if (amrex::Gpu::notInLaunchRegion())
                mfi_info.EnableTiling(amrex::IntVect(1024, 1024, 1024))
                    .SetDynamic(true);
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
            for (amrex::MFIter mfi(dof_field(lev),mfi_info); mfi.isValid(); ++mfi)
            {
                const auto& bx = mfi.tilebox();
                auto rho_arr = den(lev).array(mfi);
                auto tra_arr = dof_field(lev).array(mfi);
                amrex::FArrayBox rhotracfab;
                amrex::Array4<amrex::Real> rhotrac;

                if (PDE::multiply_rho) {
                    auto rhotrac_box = amrex::grow(bx, fvm::Godunov::nghost_state);
                    rhotracfab.resize(rhotrac_box, PDE::ndim);
                    rhotrac = rhotracfab.array();

                    amrex::ParallelFor(
                        rhotrac_box, PDE::ndim,
                        [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
                            rhotrac(i, j, k, n) =
                                rho_arr(i, j, k) * tra_arr(i, j, k, n);
                        });
                }

                amrex::FArrayBox tmpfab(amrex::grow(bx, 1), PDE::ndim * 14 + 1);

                godunov::compute_advection(lev, bx, PDE::ndim,
                                           conv_term(lev).array(mfi),
                                           (PDE::multiply_rho ? rhotrac : tra_arr),
                                           u_mac(lev).const_array(mfi),
                                           v_mac(lev).const_array(mfi),
                                           w_mac(lev).const_array(mfi),
                                           src_term(lev).const_array(mfi),
                                           dof_field.bcrec_device().data(),
                                           iconserv.data(),
                                           tmpfab.dataPtr(),
                                           geom, dt, godunov_ppm);

                amrex::Gpu::streamSynchronize();
            }
        }
    }

    PDEFields& fields;
    Field& density;
    Field& u_mac;
    Field& v_mac;
    Field& w_mac;
    amrex::Gpu::DeviceVector<int> iconserv;

    bool godunov_ppm{true};
};

} // namespace pde
} // namespace amr_wind

#endif /* ADVOP_GODUNOV_H */
