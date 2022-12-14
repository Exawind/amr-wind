#ifndef DENSITY_OPS_H
#define DENSITY_OPS_H

#include "amr-wind/equation_systems/density/density.H"

namespace amr_wind::pde {

template <typename Scheme>
struct ComputeRHSOp<Density, Scheme>
{
    explicit ComputeRHSOp(PDEFields& fields_in) : fields(fields_in) {}

    void predictor_rhs(
        const DiffusionType /*unused*/,
        const amrex::Real dt,
        bool /*mesh_mapping*/)
    {
        // Field states for diffusion and advection terms. In Godunov scheme
        // these terms only have one state.
        auto fstate = std::is_same<Scheme, fvm::Godunov>::value
                          ? FieldState::New
                          : FieldState::Old;

        const int nlevels = fields.repo.num_active_levels();
        auto& field = fields.field;
        const auto& field_old = field.state(FieldState::Old);
        const auto& conv_term = fields.conv_term.state(fstate);

        for (int lev = 0; lev < nlevels; ++lev) {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
            for (amrex::MFIter mfi(field(lev)); mfi.isValid(); ++mfi) {
                const auto& bx = mfi.tilebox();
                auto rho = field(lev).array(mfi);
                const auto rho_o = field_old(lev).const_array(mfi);
                const auto ddt_o = conv_term(lev).const_array(mfi);

                amrex::ParallelFor(
                    bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                        rho(i, j, k) = rho_o(i, j, k) + dt * ddt_o(i, j, k);
                        // Defer n+1/2 update to the time-stepping algorithm
                        // rho_nph(i, j, k) = 0.5 * (rho(i, j, k) + rho_o(i, j,
                        // k));
                    });
            }
        }
    }

    void corrector_rhs(
        const DiffusionType /*unused*/,
        const amrex::Real dt,
        bool /*mesh_mapping*/)
    {
        const int nlevels = fields.repo.num_active_levels();
        auto& field = fields.field;
        const auto& field_old = field.state(FieldState::Old);
        const auto& conv_term = fields.conv_term;
        const auto& conv_term_old = fields.conv_term.state(FieldState::Old);

        for (int lev = 0; lev < nlevels; ++lev) {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
            for (amrex::MFIter mfi(field(lev)); mfi.isValid(); ++mfi) {
                const auto& bx = mfi.tilebox();
                auto rho = field(lev).array(mfi);
                const auto rho_o = field_old(lev).const_array(mfi);
                const auto ddt = conv_term(lev).const_array(mfi);
                const auto ddt_o = conv_term_old(lev).const_array(mfi);

                amrex::ParallelFor(
                    bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                        rho(i, j, k) =
                            rho_o(i, j, k) +
                            dt * 0.5 * (ddt_o(i, j, k) + ddt(i, j, k));
                        // Defer n+1/2 update to time-stepping algorithm
                        // rho_nph(i, j, k) =
                        //     0.5 * (rho(i, j, k) + rho_o(i, j, k));
                    });
            }
        }
    }

    // data members
    PDEFields& fields;
};
} // namespace amr_wind::pde

#endif /* DENSITY_OPS_H */
