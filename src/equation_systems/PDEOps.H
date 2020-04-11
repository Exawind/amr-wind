#ifndef PDEOPS_H
#define PDEOPS_H

#include "FieldUtils.H"
#include "PDEHelpers.H"

namespace amr_wind {
namespace pde {

template<typename PDE, typename Scheme>
struct FieldRegOp
{
    FieldRegOp(FieldRepo& repo_in) : repo(repo_in) {}

    PDEFields operator()(const SimTime& time, const int probtype)
    {
        return create_fields_instance<PDE, Scheme>(time, repo, probtype);
    }

    FieldRepo& repo;
};

template<typename PDE>
struct SrcTermOpBase
{
    SrcTermOpBase(PDEFields& fields_in) : fields(fields_in) {}

    void multiply_rho(const FieldState fstate)
    {
        const auto rhostate = field_impl::phi_state(fstate);
        auto& density = fields.repo.get_field("density", rhostate);

        const int nlevels = fields.repo.num_active_levels();
        for (int lev=0; lev < nlevels; ++lev) {
            auto& src_term = fields.src_term(lev);
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (amrex::MFIter mfi(src_term, amrex::TilingIfNotGPU());
                 mfi.isValid(); ++mfi) {
                const auto& bx = mfi.tilebox();
                const auto& vf = src_term.array(mfi);
                const auto& rho = density(lev).array(mfi);

                amrex::ParallelFor(
                    bx, fields.src_term.num_comp(),
                    [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
                        vf(i, j, k, n) *= rho(i, j, k);
                    });
            }
        }
    }

    PDEFields& fields;
    amrex::Vector<std::unique_ptr<typename PDE::SrcTerm>> sources;
};

template<typename PDE>
struct SrcTermOp : SrcTermOpBase<PDE>
{
    SrcTermOp(PDEFields& fields_in) : SrcTermOpBase<PDE>(fields_in) {}

    void operator()(const FieldState fstate)
    {
        // Zero out source term
        this->fields.src_term.setVal(0.0);

        // Return early if there are no source terms to process
        if (this->sources.size() == 0) return;

        const int nlevels = this->fields.repo.num_active_levels();
        for (int lev=0; lev < nlevels; ++lev) {
            auto& src_term = this->fields.src_term(lev);
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (amrex::MFIter mfi(src_term, amrex::TilingIfNotGPU());
                 mfi.isValid(); ++mfi) {
                const auto& bx = mfi.tilebox();
                const auto& vf = src_term.array(mfi);

                for (auto& src: this->sources) {
                    (*src)(lev, mfi, bx, fstate, vf);
                }
            }
        }

        if (PDE::multiply_rho) this->multiply_rho(fstate);
    }
};

template<typename PDE, typename Scheme>
struct AdvectionOp
{};

template<typename PDE>
struct DiffusionOp
{};

template<typename PDE>
struct LinSysSolve
{};

}
}

#endif /* PDEOPS_H */
