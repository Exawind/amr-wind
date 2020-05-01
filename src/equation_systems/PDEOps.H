#ifndef PDEOPS_H
#define PDEOPS_H

#include "FieldUtils.H"
#include "PDEHelpers.H"
#include "TurbulenceModel.H"
#include "CFDSim.H"

#include "AMReX_ParmParse.H"

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
    SrcTermOpBase(PDEFields& fields_in)
    : fields(fields_in)
    , m_density(fields_in.repo.get_field("density")) {}

    void init_source_terms(const CFDSim& sim)
    {
        amrex::ParmParse pp(PDE::pde_name());
        amrex::Vector<std::string> src_terms;
        pp.queryarr("source_terms", src_terms);
        for (auto& src_name: src_terms)
        {
            sources.emplace_back(PDE::SrcTerm::create(src_name, sim));
        }
    }

    void multiply_rho(const FieldState fstate)
    {
        BL_PROFILE("amr-wind::" + PDE::pde_name() + "::multiply_rho")
        const auto rhostate = field_impl::phi_state(fstate);
        auto& density = m_density.state(rhostate);

        const int nlevels = fields.repo.num_active_levels();
        for (int lev=0; lev < nlevels; ++lev) {
            auto& src_term = fields.src_term(lev);
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
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
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
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

    PDEFields& fields;
    Field& m_density;
    amrex::Vector<std::unique_ptr<typename PDE::SrcTerm>> sources;
};

template<typename PDE>
struct SrcTermOp : SrcTermOpBase<PDE>
{
    SrcTermOp(PDEFields& fields_in) : SrcTermOpBase<PDE>(fields_in) {}
};

template<typename PDE, typename Scheme, typename = void>
struct AdvectionOp
{};

template<typename PDE, typename Scheme, typename = void>
struct DiffusionOp
{};

template<typename PDE>
struct TurbulenceOp
{
    TurbulenceOp(turbulence::TurbulenceModel& tmodel,
                 PDEFields& fields)
        : m_tmodel(tmodel), m_fields(fields)
    {}

    void operator()()
    {
        m_tmodel.update_scalar_diff(m_fields.mueff, m_fields.field.name());
    }

    turbulence::TurbulenceModel& m_tmodel;
    PDEFields& m_fields;
};

template<typename PDE, typename = void>
struct BCOp;

}
}

#endif /* PDEOPS_H */
