#ifndef VOF_OPS_H
#define VOF_OPS_H

#include "amr-wind/equation_systems/vof/vof.H"

namespace amr_wind {
namespace pde {

/** Field registration operator for Levelset equation system
 *  \ingroup levelset
 */
template <typename Scheme>
struct FieldRegOp<VOF, Scheme>
{
    explicit FieldRegOp(CFDSim& sim_in) : sim(sim_in) {}

    PDEFields operator()(const SimTime& time, const int probtype)
    {
        auto& repo = sim.repo();
        auto fields = create_fields_instance<VOF, Scheme>(time, repo, probtype);

        auto& levelset = repo.declare_cc_field("levelset", 1, 1, 1);
        auto& normal =
            repo.declare_cc_field("interface_normal", AMREX_SPACEDIM, 1, 1);
        auto& curvature = repo.declare_cc_field("interface_curvature", 1, 1, 1);

        levelset.set_default_fillpatch_bc(sim.time());
        normal.set_default_fillpatch_bc(sim.time());
        curvature.set_default_fillpatch_bc(sim.time());

        // Register fields for output/restart
        auto& iomgr = sim.io_manager();
        iomgr.register_io_var(fields.field.name());
        iomgr.register_io_var(levelset.name());
        // iomgr.register_io_var(normal.name());
        iomgr.register_io_var(curvature.name());

        return fields;
    }

    CFDSim& sim;
};

/** Right-hand side (RHS) evaluator for VOF transport equation
 */
template <typename Scheme>
struct ComputeRHSOp<VOF, Scheme>
{
    explicit ComputeRHSOp(PDEFields& fields_in) : fields(fields_in) {}

    void predictor_rhs(const DiffusionType, const amrex::Real) {}

    void corrector_rhs(const DiffusionType, const amrex::Real) {}

    // data members
    PDEFields& fields;
};
} // namespace pde
} // namespace amr_wind

#endif /* VOF_OPS_H */
