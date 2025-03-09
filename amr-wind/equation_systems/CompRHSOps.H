#ifndef COMPRHSOPS_H
#define COMPRHSOPS_H

#include "amr-wind/incflo_enums.H"
#include "amr-wind/equation_systems/PDEOps.H"
#include "amr-wind/equation_systems/SchemeTraits.H"

namespace amr_wind::pde {

/** RHS computation operator
 *  \ingroup pdeop
 *
 *  Accumulates the convective, diffusion, and source terms for predictor and
 *  corrector steps.
 */
template <typename PDE, typename Scheme>
struct ComputeRHSOp
{
    explicit ComputeRHSOp(PDEFields& fields_in)
        : fields(fields_in), density(fields_in.repo.get_field("density"))
    {}

    /** Compute right-hand side for predictor steps
     *
     *  \param difftype Indicating whether time-integration is explicit/implicit
     *  \param dt time step size
     *  \param mesh_mapping flag to indicate if mesh mapping is used
     */
    void predictor_rhs(
        const DiffusionType difftype, const amrex::Real dt, bool mesh_mapping)
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

        // for RHS evaluation velocity field should be in stretched space
        auto& field = fields.field;
        if (field.in_uniform_space() && mesh_mapping) {
            field.to_stretched_space();
        }
        auto& field_old = field.state(FieldState::Old);
        if (field_old.in_uniform_space() && mesh_mapping) {
            field_old.to_stretched_space();
        }

        auto& den_new = density.state(FieldState::New);
        auto& den_old = density.state(FieldState::Old);
        auto& src_term = fields.src_term;
        auto& diff_term = fields.diff_term.state(fstate);
        auto& conv_term = fields.conv_term.state(fstate);
        auto& mask_cell = fields.repo.get_int_field("mask_cell");
        Field const* mesh_detJ =
            mesh_mapping ? &(fields.repo.get_mesh_mapping_det_j(FieldLoc::CELL))
                         : nullptr;

        for (int lev = 0; lev < nlevels; ++lev) {
            const auto& fld_arrs = field(lev).arrays();
            const auto& fld_o_arrs = field_old(lev).const_arrays();
            const auto& rho_o_arrs = den_old(lev).const_arrays();
            const auto& rho_arrs = den_new(lev).const_arrays();
            const auto& src_arrs = src_term(lev).const_arrays();
            const auto& diff_arrs = diff_term(lev).const_arrays();
            const auto& ddt_o_arrs = conv_term(lev).const_arrays();
            const auto& imask_arrs = mask_cell(lev).const_arrays();
            const auto& detJ_arrs =
                mesh_mapping ? ((*mesh_detJ)(lev).const_arrays())
                             : amrex::MultiArray4<amrex::Real const>();

            if (PDE::multiply_rho) {
                // Remove multiplication by density as it will be added back
                // in solver
                amrex::ParallelFor(
                    field(lev), amrex::IntVect(0), PDE::ndim,
                    [=] AMREX_GPU_DEVICE(
                        int nbx, int i, int j, int k, int n) noexcept {
                        amrex::Real det_j =
                            mesh_mapping ? (detJ_arrs[nbx](i, j, k)) : 1.0;

                        fld_arrs[nbx](i, j, k, n) =
                            rho_o_arrs[nbx](i, j, k) * det_j *
                                fld_o_arrs[nbx](i, j, k, n) +
                            static_cast<amrex::Real>(imask_arrs[nbx](i, j, k)) *
                                dt *
                                (ddt_o_arrs[nbx](i, j, k, n) +
                                 det_j * src_arrs[nbx](i, j, k, n) +
                                 factor * diff_arrs[nbx](i, j, k, n));

                        fld_arrs[nbx](i, j, k, n) /= rho_arrs[nbx](i, j, k);

                        if (difftype == DiffusionType::Explicit) {
                            fld_arrs[nbx](i, j, k, n) /= det_j;
                        }
                    });
            } else {
                amrex::ParallelFor(
                    field(lev), amrex::IntVect(0), PDE::ndim,
                    [=] AMREX_GPU_DEVICE(
                        int nbx, int i, int j, int k, int n) noexcept {
                        amrex::Real det_j =
                            mesh_mapping ? (detJ_arrs[nbx](i, j, k)) : 1.0;

                        fld_arrs[nbx](i, j, k, n) =
                            det_j * fld_o_arrs[nbx](i, j, k, n) +
                            static_cast<amrex::Real>(imask_arrs[nbx](i, j, k)) *
                                dt *
                                (ddt_o_arrs[nbx](i, j, k, n) +
                                 det_j * src_arrs[nbx](i, j, k, n) +
                                 factor * diff_arrs[nbx](i, j, k, n));

                        if (difftype == DiffusionType::Explicit) {
                            fld_arrs[nbx](i, j, k, n) /= det_j;
                        }
                    });
            }
        }
        amrex::Gpu::streamSynchronize();
    }

    /** Compute right-hand side for corrector steps
     *
     *  \param difftype Indicating whether time-integration is explicit/implicit
     *  \param dt time step size
     *  \param mesh_mapping flag to indicate if mesh mapping is used
     */
    void corrector_rhs(
        const DiffusionType difftype, const amrex::Real dt, bool mesh_mapping)
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

        // for RHS evaluation velocity field should be in stretched space
        auto& field = fields.field;
        if (field.in_uniform_space() && mesh_mapping) {
            field.to_stretched_space();
        }
        auto& field_old = field.state(FieldState::Old);
        if (field_old.in_uniform_space() && mesh_mapping) {
            field_old.to_stretched_space();
        }

        auto& den_new = density.state(FieldState::New);
        auto& den_old = density.state(FieldState::Old);
        auto& src_term = fields.src_term;
        auto& diff_term = fields.diff_term;
        auto& conv_term = fields.conv_term;
        auto& diff_term_old = fields.diff_term.state(FieldState::Old);
        auto& conv_term_old = fields.conv_term.state(FieldState::Old);
        auto& mask_cell = fields.repo.get_int_field("mask_cell");
        Field const* mesh_detJ =
            mesh_mapping ? &(fields.repo.get_mesh_mapping_det_j(FieldLoc::CELL))
                         : nullptr;

        for (int lev = 0; lev < nlevels; ++lev) {
            const auto& fld_arrs = field(lev).arrays();
            const auto& fld_o_arrs = field_old(lev).const_arrays();
            const auto& rho_o_arrs = den_old(lev).const_arrays();
            const auto& rho_arrs = den_new(lev).const_arrays();
            const auto& src_arrs = src_term(lev).const_arrays();
            const auto& diff_arrs = diff_term(lev).const_arrays();
            const auto& ddt_arrs = conv_term(lev).const_arrays();
            const auto& diff_o_arrs = diff_term_old(lev).const_arrays();
            const auto& ddt_o_arrs = conv_term_old(lev).const_arrays();
            const auto& imask_arrs = mask_cell(lev).const_arrays();
            const auto& detJ_arrs =
                mesh_mapping ? ((*mesh_detJ)(lev).const_arrays())
                             : amrex::MultiArray4<amrex::Real const>();

            if (PDE::multiply_rho) {
                // Remove multiplication by density as it will be added back
                // in solver
                amrex::ParallelFor(
                    field(lev), amrex::IntVect(0), PDE::ndim,
                    [=] AMREX_GPU_DEVICE(
                        int nbx, int i, int j, int k, int n) noexcept {
                        amrex::Real det_j =
                            mesh_mapping ? (detJ_arrs[nbx](i, j, k)) : 1.0;

                        fld_arrs[nbx](i, j, k, n) =
                            rho_o_arrs[nbx](i, j, k) * det_j *
                                fld_o_arrs[nbx](i, j, k, n) +
                            static_cast<amrex::Real>(imask_arrs[nbx](i, j, k)) *
                                dt *
                                (0.5 * (ddt_o_arrs[nbx](i, j, k, n) +
                                        ddt_arrs[nbx](i, j, k, n)) +
                                 ofac * diff_o_arrs[nbx](i, j, k, n) +
                                 nfac * diff_arrs[nbx](i, j, k, n) +
                                 det_j * src_arrs[nbx](i, j, k, n));

                        fld_arrs[nbx](i, j, k, n) /= rho_arrs[nbx](i, j, k);

                        if (difftype == DiffusionType::Explicit) {
                            fld_arrs[nbx](i, j, k, n) /= det_j;
                        }
                    });
            } else {
                amrex::ParallelFor(
                    field(lev), amrex::IntVect(0), PDE::ndim,
                    [=] AMREX_GPU_DEVICE(
                        int nbx, int i, int j, int k, int n) noexcept {
                        amrex::Real det_j =
                            mesh_mapping ? (detJ_arrs[nbx](i, j, k)) : 1.0;

                        fld_arrs[nbx](i, j, k, n) =
                            det_j * fld_o_arrs[nbx](i, j, k, n) +
                            static_cast<amrex::Real>(imask_arrs[nbx](i, j, k)) *
                                dt *
                                (0.5 * (ddt_o_arrs[nbx](i, j, k, n) +
                                        ddt_arrs[nbx](i, j, k, n)) +
                                 ofac * diff_o_arrs[nbx](i, j, k, n) +
                                 nfac * diff_arrs[nbx](i, j, k, n) +
                                 det_j * src_arrs[nbx](i, j, k, n));

                        if (difftype == DiffusionType::Explicit) {
                            fld_arrs[nbx](i, j, k, n) /= det_j;
                        }
                    });
            }
        }
        amrex::Gpu::streamSynchronize();
    }

    // data members
    PDEFields& fields;
    Field& density;
};

} // namespace amr_wind::pde

#endif /* COMPRHSOPS_H */
