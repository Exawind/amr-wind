#include "DiffusionOps.H"
#include "diffusion.H"
#include "console_io.H"

#include "AMReX_MLTensorOp.H"

namespace amr_wind {
namespace pde {

template <typename LinOp>
DiffSolverIface<LinOp>::DiffSolverIface(
    PDEFields& fields, const std::string& prefix)
    : m_pdefields(fields)
    , m_density(fields.repo.get_field("density"))
    , m_options(prefix)
{
    amrex::LPInfo isolve;
    amrex::LPInfo iapply;

    isolve.setMaxCoarseningLevel(m_options.max_coarsen_level);
    iapply.setMaxCoarseningLevel(0);

    const auto& mesh = m_pdefields.repo.mesh();
    m_solver.reset(new LinOp(mesh.Geom(0, mesh.finestLevel()),
                             mesh.boxArray(0, mesh.finestLevel()),
                             mesh.DistributionMap(0, mesh.finestLevel()),
                             isolve));
    m_applier.reset(new LinOp(mesh.Geom(0, mesh.finestLevel()),
                             mesh.boxArray(0, mesh.finestLevel()),
                             mesh.DistributionMap(0, mesh.finestLevel()),
                             iapply));

    m_solver->setMaxOrder(m_options.max_order);
    m_applier->setMaxOrder(m_options.max_order);

    // It is the sub-classes responsibility to set the linear solver BC for the
    // operators.
}

template <typename LinOp>
void DiffSolverIface<LinOp>::setup_operator(
    LinOp& linop, const amrex::Real alpha, const amrex::Real beta, const FieldState fstate)
{
    BL_PROFILE("amr-wind::setup_operator")
    auto& repo = m_pdefields.repo;
    const int nlevels = repo.num_active_levels();
    auto& density = m_density.state(fstate);

    linop.setScalars(alpha, beta);
    for (int lev = 0; lev < nlevels; ++lev) {
        linop.setACoeffs(lev, density(lev));
        linop.setLevelBC(lev, &m_pdefields.field(lev));
    }
    set_bcoeffs(linop);
}

template<typename LinOp>
void DiffSolverIface<LinOp>::setup_solver(amrex::MLMG& mlmg)
{
    BL_PROFILE("amr-wind::setup_solver")
    auto& opts = m_options;

    if (opts.bottom_solver_type == "smoother") {
        mlmg.setBottomSolver(amrex::MLMG::BottomSolver::smoother);
    } else if (opts.bottom_solver_type == "hypre") {
        mlmg.setBottomSolver(amrex::MLMG::BottomSolver::hypre);
    }

    mlmg.setMaxIter(opts.max_iter);
    mlmg.setMaxFmgIter(opts.fmg_max_iter);
    mlmg.setCGMaxIter(opts.cg_max_iter);
    mlmg.setVerbose(opts.verbose);
    mlmg.setCGVerbose(opts.cg_verbose);
}


template<typename LinOp>
void DiffSolverIface<LinOp>::linsys_solve(const amrex::Real dt)
{
    FieldState fstate = FieldState::New;
    this->setup_operator(*this->m_solver, 1.0, dt, fstate);

    auto& repo = this->m_pdefields.repo;
    auto& field = this->m_pdefields.field;
    auto& density = m_density.state(fstate);
    const int nlevels = repo.num_active_levels();
    const int ndim = field.num_comp();
    auto rhs_ptr = repo.create_scratch_field("rhs", field.num_comp(), 0);

    // Always multiply with rho since there is no diffusion term for density
    for (int lev = 0; lev < nlevels; ++lev) {
        auto& rhs = (*rhs_ptr)(lev);

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        for (amrex::MFIter mfi(rhs, amrex::TilingIfNotGPU()); mfi.isValid();
             ++mfi) {
            const auto& bx = mfi.tilebox();
            const auto& rhs_a = rhs.array(mfi);
            const auto& fld = field(lev).const_array(mfi);
            const auto& rho = density(lev).const_array(mfi);

            amrex::ParallelFor(
                bx, ndim,
                [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
                    rhs_a(i, j, k, n) = rho(i, j, k) * fld(i, j, k, n);
                });
        }
    }

    amrex::MLMG mlmg(*this->m_solver);
    this->setup_solver(mlmg);

    mlmg.solve(
        field.vec_ptrs(), rhs_ptr->vec_const_ptrs(), this->m_options.rel_tol,
        this->m_options.abs_tol);

    io::print_mlmg_info(field.name() + "_solve", mlmg);
}

template class DiffSolverIface<amrex::MLABecLaplacian>;
template class DiffSolverIface<amrex::MLTensorOp>;


} // namespace pde
} // namespace amr_wind
