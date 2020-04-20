#ifndef PDE_H
#define PDE_H

#include <string>

#include "incflo_enums.H"
#include "Factory.H"
#include "PDEHelpers.H"
#include "PDEOps.H"
#include "CompRHSOps.H"
#include "DiffusionOps.H"

namespace amr_wind {
namespace pde {

class PDEBase : public Factory<PDEBase, const SimTime&, FieldRepo&, const int>
{
public:
    virtual ~PDEBase() = default;

    virtual PDEFields& fields() = 0;

    virtual void initialize() = 0;

    virtual void post_regrid_actions() = 0;

    virtual void compute_source_term(const FieldState fstate) = 0;

    virtual void compute_nueff(const FieldState fstate) = 0;

    virtual void compute_diffusion_term(const FieldState fstate) = 0;

    virtual void compute_advection_term(const FieldState fstate) = 0;

    virtual void compute_predictor_rhs(const DiffusionType difftype) = 0;

    virtual void compute_corrector_rhs(const DiffusionType difftype) = 0;

    virtual void solve(const amrex::Real dt) = 0;

    static std::string base_identifier()
    {
        return "PDESystem";
    }
};

template<typename PDE, typename Scheme>
class PDESystem : public PDEBase::Register<PDESystem<PDE, Scheme>>
{
public:
    using PDEType = PDE;
    using SchemeType = Scheme;

    static std::string identifier()
    {
        return PDE::pde_name() + "-"  + Scheme::scheme_name();
    }

    PDESystem(const SimTime& time, FieldRepo& repo, const int probtype)
        : m_time(time)
        , m_repo(repo)
        , m_fields(FieldRegOp<PDE, Scheme>(repo)(time, probtype))
        , m_src_op(m_fields)
        , m_adv_op(m_fields)
        , m_rhs_op(m_fields)
    {}

    void initialize() override
    {
        if (PDE::has_diffusion)
            m_diff_op.reset(new DiffusionOp<PDE, Scheme>(m_fields));
    }

    void post_regrid_actions() override
    {
        if (PDE::has_diffusion)
            m_diff_op.reset(new DiffusionOp<PDE, Scheme>(m_fields));
    }

    PDEFields& fields() override { return m_fields; }

    void compute_source_term(const FieldState fstate) override
    {
        m_src_op(fstate);
    }

    void compute_nueff(const FieldState ) override
    {
        // TBD
    }

    void compute_diffusion_term(const FieldState fstate) override
    {
        if (PDE::has_diffusion)
            m_diff_op->compute_diff_term(fstate);
    }

    void compute_advection_term(const FieldState fstate) override
    {
        m_adv_op(fstate, m_time.deltaT());
    }

    virtual void compute_predictor_rhs(const DiffusionType difftype) override
    {
        m_rhs_op.predictor_rhs(difftype, m_time.deltaT());
    }

    virtual void compute_corrector_rhs(const DiffusionType difftype) override
    {
        m_rhs_op.corrector_rhs(difftype, m_time.deltaT());
    }

    void solve(const amrex::Real dt) override
    {
        if (PDE::has_diffusion)
            m_diff_op->linsys_solve(dt);
    }

protected:
    const SimTime& m_time;

    FieldRepo& m_repo;

    PDEFields m_fields;

    SrcTermOp<PDE> m_src_op;
    AdvectionOp<PDE, Scheme> m_adv_op;
    ComputeRHSOp<PDE, Scheme> m_rhs_op;
    std::unique_ptr<DiffusionOp<PDE, Scheme>> m_diff_op;
};

} // namespace pde
} // namespace amr_wind

#endif /* PDE_H */
