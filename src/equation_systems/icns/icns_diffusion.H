#ifndef ICNS_DIFFUSION_H
#define ICNS_DIFFUSION_H

#include "PDETraits.H"
#include "PDEOps.H"
#include "PDEHelpers.H"
#include "DiffusionOps.H"
#include "icns/icns.H"

namespace amr_wind {
namespace pde {

template <typename Scheme>
struct DiffusionOp<ICNS, Scheme>
    : public DiffSolverIface<typename ICNS::MLDiffOp>
{
    static_assert(
        ICNS::ndim == AMREX_SPACEDIM,
        "DiffusionOp invoked for non-scalar PDE type");
    static_assert(
        std::is_same<typename ICNS::MLDiffOp, amrex::MLTensorOp>::value,
        "Invalid linear operator for scalar diffusion operator");

    DiffusionOp(PDEFields& fields)
        : DiffSolverIface<typename ICNS::MLDiffOp>(fields)
    {
        this->m_solver->setDomainBC(
            diffusion::get_diffuse_tensor_bc(
                this->m_pdefields.field, amrex::Orientation::low),
            diffusion::get_diffuse_tensor_bc(
                this->m_pdefields.field, amrex::Orientation::high));
        this->m_applier->setDomainBC(
            diffusion::get_diffuse_tensor_bc(
                this->m_pdefields.field, amrex::Orientation::low),
            diffusion::get_diffuse_tensor_bc(
                this->m_pdefields.field, amrex::Orientation::high));
    }

    void compute_diff_term(const FieldState fstate)
    {
        this->setup_operator(*this->m_applier, 0.0, -1.0, fstate);

        auto tau_state = std::is_same<Scheme, fvm::Godunov>::value
                             ? FieldState::New
                             : fstate;
        auto& divtau = this->m_pdefields.diff_term.state(tau_state);

        amrex::MLMG mlmg(*this->m_applier);
        mlmg.apply(divtau.vec_ptrs(), this->m_pdefields.field.vec_ptrs());

        const auto& repo = this->m_pdefields.repo;
        const int nlevels = repo.num_active_levels();
        const auto& density = repo.get_field("density", fstate);

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        for (int lev = 0; lev < nlevels; ++lev) {
            for (amrex::MFIter mfi(divtau(lev), amrex::TilingIfNotGPU());
                 mfi.isValid(); ++mfi) {
                const auto& bx = mfi.tilebox();
                const auto& divtau_arr = divtau(lev).array(mfi);
                const auto& rho_arr = density(lev).const_array(mfi);

                amrex::ParallelFor(
                    bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                        amrex::Real rhoinv = 1.0 / rho_arr(i, j, k);
                        divtau_arr(i, j, k, 0) *= rhoinv;
                        divtau_arr(i, j, k, 1) *= rhoinv;
                        divtau_arr(i, j, k, 2) *= rhoinv;
                    });
            }
        }
    }
};

} // namespace pde
} // namespace amr_wind

#endif /* ICNS_DIFFUSION_H */
