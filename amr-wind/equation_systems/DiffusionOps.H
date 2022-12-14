#ifndef DIFFUSIONOPS_H
#define DIFFUSIONOPS_H

#include "amr-wind/core/MLMGOptions.H"
#include "amr-wind/equation_systems/PDETraits.H"
#include "amr-wind/equation_systems/SchemeTraits.H"
#include "amr-wind/equation_systems/PDEOps.H"
#include "amr-wind/equation_systems/PDEHelpers.H"
#include "amr-wind/diffusion/diffusion.H"

#include "AMReX_MLABecLaplacian.H"
#include "AMReX_MLTensorOp.H"
#include "AMReX_MLMG.H"
#include "AMReX_AmrCore.H"
#include "AMReX_MultiFabUtil.H"

namespace amr_wind::pde {

/** Base class for all diffusion operators
 *  \ingroup pdeop
 *
 *  This class provides the common operations for an implicit solution of a
 *  convection-diffusion equation within AMR-Wind.
 *
 *  \tparam LinOp The linear operator (see [AMREeX
 * docs](https://amrex-codes.github.io/amrex/docs_html/LinearSolvers.html))
 */
template <typename LinOp>
class DiffSolverIface
{
public:
    DiffSolverIface(
        PDEFields& fields,
        const bool has_overset,
        const bool mesh_mapping,
        const std::string& prefix = "diffusion");

    virtual ~DiffSolverIface() = default;

    /** Implicit solve and update of a linear system
     *
     *  \param dt timestep size
     */
    virtual void linsys_solve(const amrex::Real dt);

    virtual void linsys_solve_impl();

    virtual void set_acoeffs(LinOp& linop, const FieldState fstate);

    template <typename L>
    void set_bcoeffs(
        L& linop,
        typename std::enable_if<
            std::is_same<L, amrex::MLTensorOp>::value>::type* /*unused*/
        = nullptr)
    {
        const int nlevels = m_pdefields.repo.num_active_levels();
        const auto& viscosity = m_pdefields.mueff;
        const auto& geom = m_pdefields.repo.mesh().Geom();

        for (int lev = 0; lev < nlevels; ++lev) {
            auto b = diffusion::average_velocity_eta_to_faces(
                geom[lev], viscosity(lev));
            if (m_mesh_mapping) {
                diffusion::viscosity_to_uniform_space(b, m_pdefields.repo, lev);
            }
            linop.setShearViscosity(lev, amrex::GetArrOfConstPtrs(b));
        }
    }

    template <typename L>
    void set_bcoeffs(
        L& linop,
        typename std::enable_if<
            std::is_same<L, amrex::MLABecLaplacian>::value>::type* /*unused*/
        = nullptr)
    {
        const int nlevels = m_pdefields.repo.num_active_levels();
        const auto& viscosity = m_pdefields.mueff;
        const auto& geom = m_pdefields.repo.mesh().Geom();

        for (int lev = 0; lev < nlevels; ++lev) {
            auto b = diffusion::average_velocity_eta_to_faces(
                geom[lev], viscosity(lev));
            if (m_mesh_mapping) {
                diffusion::viscosity_to_uniform_space(b, m_pdefields.repo, lev);
            }
            linop.setBCoeffs(lev, amrex::GetArrOfConstPtrs(b));
        }
    }

protected:
    //! Sets up the linear operator (e.g., setup BCs, etc.)
    virtual void setup_operator(
        LinOp& linop,
        const amrex::Real alpha,
        const amrex::Real beta,
        const FieldState fstate);

    virtual void setup_solver(amrex::MLMG& mlmg);

    PDEFields& m_pdefields;
    Field& m_density;

    MLMGOptions m_options;

    bool m_mesh_mapping{false};

    std::unique_ptr<LinOp> m_solver;
    std::unique_ptr<LinOp> m_applier;
};

/** Diffusion operator for scalar transport equations
 *  \ingroup pdeop
 */
template <typename PDE, typename Scheme>
struct DiffusionOp<
    PDE,
    Scheme,
    typename std::enable_if<std::is_base_of<ScalarTransport, PDE>::value>::type>
    : public DiffSolverIface<typename PDE::MLDiffOp>
{
    static_assert(
        PDE::ndim == 1, "DiffusionOp invoked for non-scalar PDE type");
    static_assert(
        std::is_same<typename PDE::MLDiffOp, amrex::MLABecLaplacian>::value,
        "Invalid linear operator for scalar diffusion operator");

    DiffusionOp(
        PDEFields& fields, const bool has_overset, const bool mesh_mapping)
        : DiffSolverIface<typename PDE::MLDiffOp>(
              fields, has_overset, mesh_mapping)
    {
        this->m_solver->setDomainBC(
            diffusion::get_diffuse_scalar_bc(
                this->m_pdefields.field, amrex::Orientation::low),
            diffusion::get_diffuse_scalar_bc(
                this->m_pdefields.field, amrex::Orientation::high));
        this->m_applier->setDomainBC(
            diffusion::get_diffuse_scalar_bc(
                this->m_pdefields.field, amrex::Orientation::low),
            diffusion::get_diffuse_scalar_bc(
                this->m_pdefields.field, amrex::Orientation::high));
    }

    //! Computes the diffusion term that goes in the RHS
    void compute_diff_term(const FieldState fstate)
    {
        this->setup_operator(*this->m_applier, 0.0, -1.0, fstate);

        auto tau_state = std::is_same<Scheme, fvm::Godunov>::value
                             ? FieldState::New
                             : fstate;
        amrex::MLMG mlmg(*this->m_applier);
        mlmg.apply(
            this->m_pdefields.diff_term.state(tau_state).vec_ptrs(),
            this->m_pdefields.field.vec_ptrs());
    }
};

} // namespace amr_wind::pde

#endif /* DIFFUSIONOPS_H */
