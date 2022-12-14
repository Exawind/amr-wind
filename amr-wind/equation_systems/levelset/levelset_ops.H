#ifndef LEVELSET_OPS_H
#define LEVELSET_OPS_H

#include "amr-wind/equation_systems/levelset/levelset.H"
#include "amr-wind/physics/multiphase/MultiPhase.H"

namespace amr_wind::pde {

/** Field registration operator for Levelset equation system
 *  \ingroup levelset
 */
template <typename Scheme>
struct FieldRegOp<Levelset, Scheme>
{
    explicit FieldRegOp(CFDSim& sim_in) : sim(sim_in) {}

    PDEFields operator()(const SimTime& time)
    {
        auto& repo = sim.repo();
        auto fields = create_fields_instance<Levelset, Scheme>(time, repo);

        auto& normal =
            repo.declare_cc_field("interface_normal", AMREX_SPACEDIM, 1, 1);
        auto& curvature = repo.declare_cc_field("interface_curvature", 1, 1, 1);

        normal.template register_fill_patch_op<FieldFillPatchOps<FieldBCNoOp>>(
            repo.mesh(), time);
        curvature
            .template register_fill_patch_op<FieldFillPatchOps<FieldBCNoOp>>(
                repo.mesh(), time);

        // Register fields for output/restart
        auto& iomgr = sim.io_manager();
        iomgr.register_io_var(fields.field.name());
        iomgr.register_io_var(normal.name());
        iomgr.register_io_var(curvature.name());

        return fields;
    }

    CFDSim& sim;
};

/** Right-hand side (RHS) evaluator for Levelset transport equation
 *  \ingroup levelset
 */
template <typename Scheme>
struct ComputeRHSOp<Levelset, Scheme>
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
                auto phi = field(lev).array(mfi);
                const auto phi_o = field_old(lev).const_array(mfi);
                const auto ddt_o = conv_term(lev).const_array(mfi);

                amrex::ParallelFor(
                    bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                        phi(i, j, k) = phi_o(i, j, k) + dt * ddt_o(i, j, k);
                        // Defer n+1/2 update to the time-stepping algorithm
                        // phi_nph(i, j, k) = 0.5 * (phi(i, j, k) + phi_o(i, j,
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
                auto phi = field(lev).array(mfi);
                const auto phi_o = field_old(lev).const_array(mfi);
                const auto ddt = conv_term(lev).const_array(mfi);
                const auto ddt_o = conv_term_old(lev).const_array(mfi);

                amrex::ParallelFor(
                    bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                        phi(i, j, k) =
                            phi_o(i, j, k) +
                            dt * 0.5 * (ddt_o(i, j, k) + ddt(i, j, k));
                        // Defer n+1/2 update to time-stepping algorithm
                        // phi_nph(i, j, k) =
                        //     0.5 * (phi(i, j, k) + phi_o(i, j, k));
                    });
            }
        }
    }

    // data members
    PDEFields& fields;
};

/** Implementation of post-processing actions for the levelset equation
 *  \ingroup levelset
 */
template <>
struct PostSolveOp<Levelset>
{
    PostSolveOp(CFDSim& sim, PDEFields& fields) : m_sim(sim), m_fields(fields)
    {}

    void operator()(const amrex::Real /*unused*/)
    {
        m_fields.field.fillpatch(m_sim.time().current_time());
        auto& multiphase = m_sim.physics_manager().get<MultiPhase>();
        multiphase.set_density_via_levelset();
    }

    CFDSim& m_sim;
    PDEFields& m_fields;
};

} // namespace amr_wind::pde

#endif /* LEVELSET_OPS_H */
