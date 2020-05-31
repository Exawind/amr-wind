#ifndef COMPRHSOPS_H
#define COMPRHSOPS_H

#include "amr-wind/incflo_enums.H"
#include "amr-wind/equation_systems/PDEOps.H"
#include "amr-wind/equation_systems/SchemeTraits.H"

namespace amr_wind {
namespace pde {

template<typename PDE, typename Scheme>
struct ComputeRHSOp
{
    ComputeRHSOp(PDEFields& fields_in)
    : fields(fields_in)
    , density(fields_in.repo.get_field("density"))
    {}

    void predictor_rhs(const DiffusionType difftype, const amrex::Real dt)
    {
        amrex::Real factor = 0.0;
        switch (difftype) {
        case DiffusionType::Explicit:
            factor = 1.0;
            break;

        case DiffusionType::Crank_Nicolson:
            factor = 0.5;
            break;

        case DiffusionType::Implicit:
            factor = 0.0;
            break;

        default:
            amrex::Abort("Invalid diffusion type");
        }

        // Field states for diffusion and advection terms. In Godunov scheme
        // these terms only have one state.
        auto fstate = std::is_same<Scheme, fvm::Godunov>::value
                          ? FieldState::New
                          : FieldState::Old;

        const int nlevels = fields.repo.num_active_levels();
        auto& field = fields.field;
        auto& field_old = field.state(FieldState::Old);
        auto& den_new = density.state(FieldState::New);
        auto& den_old = density.state(FieldState::Old);
        auto& src_term = fields.src_term;
        auto& diff_term = fields.diff_term.state(fstate);
        auto& conv_term = fields.conv_term.state(fstate);

        for (int lev = 0; lev < nlevels; ++lev) {
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
            for (amrex::MFIter mfi(field(lev)); mfi.isValid(); ++mfi) {
                const auto& bx = mfi.tilebox();
                auto fld = field(lev).array(mfi);
                const auto fld_o = field_old(lev).const_array(mfi);
                const auto rho_o = den_old(lev).const_array(mfi);
                const auto rho = den_new(lev).const_array(mfi);
                const auto src = src_term(lev).const_array(mfi);
                const auto diff = diff_term(lev).const_array(mfi);
                const auto ddt_o = conv_term(lev).const_array(mfi);

                if (PDE::multiply_rho) {
                    // Remove multiplication by density as it will be added back
                    // in solver
                    amrex::ParallelFor(
                        bx, PDE::ndim,
                        [=] AMREX_GPU_DEVICE(
                            int i, int j, int k, int n) noexcept {
                            fld(i, j, k, n) =
                                rho_o(i, j, k) * fld_o(i, j, k, n) +
                                dt * (ddt_o(i, j, k, n) + src(i, j, k, n) +
                                      factor * diff(i, j, k, n));

                            fld(i, j, k, n) /= rho(i, j, k);
                        });
                } else {
                    amrex::ParallelFor(
                        bx, PDE::ndim,
                        [=] AMREX_GPU_DEVICE(
                            int i, int j, int k, int n) noexcept {
                            fld(i, j, k, n) =
                                fld_o(i, j, k, n) +
                                dt * (ddt_o(i, j, k, n) + src(i, j, k, n) +
                                      factor * diff(i, j, k, n));
                        });
                }
            }
        }
    }

    void corrector_rhs(const DiffusionType difftype, const amrex::Real dt)
    {
        amrex::Real ofac = 0.0;
        amrex::Real nfac = 0.0;
        switch (difftype) {
        case DiffusionType::Explicit:
            ofac = 0.5;
            nfac = 0.5;
            break;

        case DiffusionType::Crank_Nicolson:
            ofac = 0.5;
            nfac = 0.0;
            break;

        case DiffusionType::Implicit:
            ofac = 0.0;
            nfac = 0.0;
            break;

        default:
            amrex::Abort("Invalid diffusion type");
        }

        const int nlevels = fields.repo.num_active_levels();
        auto& field = fields.field;
        auto& field_old = field.state(FieldState::Old);
        auto& den_new = density.state(FieldState::New);
        auto& den_old = density.state(FieldState::Old);
        auto& src_term = fields.src_term;
        auto& diff_term = fields.diff_term;
        auto& conv_term = fields.conv_term;
        auto& diff_term_old = fields.diff_term.state(FieldState::Old);
        auto& conv_term_old = fields.conv_term.state(FieldState::Old);

        for (int lev = 0; lev < nlevels; ++lev) {
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
            for (amrex::MFIter mfi(field(lev)); mfi.isValid(); ++mfi) {
                const auto& bx = mfi.tilebox();
                auto fld = field(lev).array(mfi);
                const auto fld_o = field_old(lev).const_array(mfi);
                const auto rho_o = den_old(lev).const_array(mfi);
                const auto rho = den_new(lev).const_array(mfi);
                const auto src = src_term(lev).const_array(mfi);
                const auto diff = diff_term(lev).const_array(mfi);
                const auto ddt = conv_term(lev).const_array(mfi);
                const auto diff_o = diff_term_old(lev).const_array(mfi);
                const auto ddt_o = conv_term_old(lev).const_array(mfi);

                if (PDE::multiply_rho) {
                    // Remove multiplication by density as it will be added back
                    // in solver
                    amrex::ParallelFor(
                        bx, PDE::ndim,
                        [=] AMREX_GPU_DEVICE(
                            int i, int j, int k, int n) noexcept {
                            fld(i, j, k, n) =
                                rho_o(i, j, k) * fld_o(i, j, k, n) +
                                dt *
                                    (0.5 *
                                         (ddt_o(i, j, k, n) + ddt(i, j, k, n)) +
                                     ofac * diff_o(i, j, k, n) +
                                     nfac * diff(i, j, k, n) + src(i, j, k, n));

                            fld(i, j, k, n) /= rho(i, j, k);
                        });
                } else {
                    amrex::ParallelFor(
                        bx, PDE::ndim,
                        [=] AMREX_GPU_DEVICE(
                            int i, int j, int k, int n) noexcept {
                            fld(i, j, k, n) = fld_o(i, j, k, n) +
                                              dt * (0.5 * (ddt_o(i, j, k, n) +
                                                           ddt(i, j, k, n)) +
                                                    ofac * diff_o(i, j, k, n) +
                                                    nfac * diff(i, j, k, n) +
                                                    src(i, j, k, n));
                        });
                }
            }
        }
    }

    // data members
    PDEFields& fields;
    Field& density;
};

}
}

#endif /* COMPRHSOPS_H */
