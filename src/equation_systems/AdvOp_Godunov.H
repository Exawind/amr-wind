#ifndef ADVOP_GODUNOV_H
#define ADVOP_GODUNOV_H

#include <type_traits>

#include "Godunov.H"
#include "SchemeTraits.H"
#include "PDETraits.H"
#include "PDEOps.H"

#include "AMReX_Gpu.H"
#include "AMReX_ParmParse.H"

namespace amr_wind {
namespace pde {

template <typename PDE>
struct AdvectionOp<
    PDE,
    fvm::Godunov,
    typename std::enable_if<std::is_base_of<ScalarTransport, PDE>::value>::type>
{
    AdvectionOp(PDEFields& fields_in) : fields(fields_in)
    {
        amrex::ParmParse pp("incflo");
        pp.query("godunov_ppm", godunov_ppm);

        // TODO: Need iconserv flag to be adjusted???
        iconserv.resize(PDE::ndim, 1);
    }

    void operator()(const FieldState fstate, const amrex::Real dt)
    {
        static_assert(PDE::ndim == 1, "Invalid number of components for scalar");
        auto& repo = fields.repo;
        auto& geom = repo.mesh().Geom();
        auto& u_mac = repo.get_field("u_mac");
        auto& v_mac = repo.get_field("v_mac");
        auto& w_mac = repo.get_field("w_mac");

        auto& src_term = fields.src_term;
        auto& conv_term = fields.conv_term;
        auto& dof_field = fields.field.state(fstate);
        auto& den = repo.get_field("density", fstate);

        for (int lev = 0; lev < repo.num_active_levels(); ++lev) {
            amrex::MFItInfo mfi_info;
            if (amrex::Gpu::notInLaunchRegion())
                mfi_info.EnableTiling(amrex::IntVect(1024, 1024, 1024))
                    .SetDynamic(true);
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (amrex::MFIter mfi(den(lev),mfi_info); mfi.isValid(); ++mfi)
            {
                const auto& bx = mfi.tilebox();
                auto rho_arr = den(lev).array(mfi);
                auto tra_arr = dof_field(lev).array(mfi);
                amrex::Array4<amrex::Real> rhotrac;

                if (PDE::multiply_rho) {
                    auto rhotrac_box = amrex::grow(bx, fvm::Godunov::nghost_state);
                    amrex::FArrayBox rhotracfab(rhotrac_box, PDE::ndim);
                    rhotrac = rhotracfab.array();

                    amrex::ParallelFor(
                        rhotrac_box,
                        [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                            rhotrac(i, j, k) =
                                rho_arr(i, j, k) * tra_arr(i, j, k);
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
    amrex::Gpu::DeviceVector<int> iconserv;

    bool godunov_ppm{true};
};

} // namespace pde
} // namespace amr_wind

#endif /* ADVOP_GODUNOV_H */
