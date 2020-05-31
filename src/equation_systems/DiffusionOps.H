#ifndef DIFFUSIONOPS_H
#define DIFFUSIONOPS_H

#include "amr-wind/core/MLMGOptions.H"
#include "amr-wind/equation_systems/PDETraits.H"
#include "amr-wind/equation_systems/SchemeTraits.H"
#include "amr-wind/equation_systems/PDEOps.H"
#include "amr-wind/equation_systems/PDEHelpers.H"
#include "amr-wind/diffusion/diffusion.H"

#include "AMReX_MLABecLaplacian.H"
#include "AMReX_MLMG.H"
#include "AMReX_AmrCore.H"

namespace amr_wind {
namespace pde {

template <typename LinOp>
class DiffSolverIface
{
public:
    DiffSolverIface(PDEFields& fields, const std::string& prefix = "diffusion");

    virtual ~DiffSolverIface() = default;

    virtual void linsys_solve(const amrex::Real dt);

protected:
    virtual void setup_operator(
        LinOp& linop,
        const amrex::Real alpha,
        const amrex::Real beta,
        const FieldState fstate);

    virtual void setup_solver(amrex::MLMG& mlmg);

    template <typename L>
    void set_bcoeffs(
        L& linop,
        typename std::enable_if<
            std::is_same<L, amrex::MLTensorOp>::value>::type* = 0)
    {
        const int nlevels = m_pdefields.repo.num_active_levels();
        auto& viscosity = m_pdefields.mueff;
        auto& geom = m_pdefields.repo.mesh().Geom();
        for (int lev = 0; lev < nlevels; ++lev) {
            auto b = diffusion::average_velocity_eta_to_faces(
                geom[lev], viscosity(lev));
            linop.setShearViscosity(lev, amrex::GetArrOfConstPtrs(b));
        }
    }

    template <typename L>
    void set_bcoeffs(
        L& linop,
        typename std::enable_if<
            std::is_same<L, amrex::MLABecLaplacian>::value>::type* = 0)
    {
        const int nlevels = m_pdefields.repo.num_active_levels();
        auto& viscosity = m_pdefields.mueff;
        auto& geom = m_pdefields.repo.mesh().Geom();
        for (int lev = 0; lev < nlevels; ++lev) {
            auto b = diffusion::average_velocity_eta_to_faces(
                geom[lev], viscosity(lev));
            linop.setBCoeffs(lev, amrex::GetArrOfConstPtrs(b));
        }
    }

    PDEFields& m_pdefields;
    Field& m_density;
    
    MLMGOptions m_options;

    std::unique_ptr<LinOp> m_solver;
    std::unique_ptr<LinOp> m_applier;
};

template <typename PDE, typename Scheme>
struct DiffusionOp<
    PDE,
    Scheme,
    typename std::enable_if<std::is_base_of<ScalarTransport, PDE>::value>::type>
    : public DiffSolverIface<typename PDE::MLDiffOp>
{
    static_assert(PDE::ndim == 1, "DiffusionOp invoked for non-scalar PDE type");
    static_assert(std::is_same<typename PDE::MLDiffOp, amrex::MLABecLaplacian>::value,
                  "Invalid linear operator for scalar diffusion operator");

    DiffusionOp(PDEFields& fields)
        : DiffSolverIface<typename PDE::MLDiffOp>(fields)
    {
        this->m_solver->setDomainBC(
            diffusion::get_diffuse_scalar_bc(this->m_pdefields.field, amrex::Orientation::low),
            diffusion::get_diffuse_scalar_bc(this->m_pdefields.field, amrex::Orientation::high));
        this->m_applier->setDomainBC(
            diffusion::get_diffuse_scalar_bc(this->m_pdefields.field, amrex::Orientation::low),
            diffusion::get_diffuse_scalar_bc(this->m_pdefields.field, amrex::Orientation::high));
    }

    void compute_diff_term(const FieldState fstate)
    {
        this->setup_operator(*this->m_applier, 0.0, -1.0, fstate);

        auto tau_state = std::is_same<Scheme, fvm::Godunov>::value
            ? FieldState::New : fstate;
        amrex::MLMG mlmg(*this->m_applier);
        mlmg.apply(this->m_pdefields.diff_term.state(tau_state).vec_ptrs(),
                   this->m_pdefields.field.vec_ptrs());
    }
};

} // namespace pde
} // namespace amr_wind

#endif /* DIFFUSIONOPS_H */
