#ifndef PDE_H
#define PDE_H

#include <string>

#include "Factory.H"
#include "PDEHelpers.H"
#include "PDEOps.H"

namespace amr_wind {
namespace pde {

class PDEBase : public Factory<PDEBase, const SimTime&, FieldRepo&, const int>
{
public:
    virtual ~PDEBase() = default;

    virtual PDEFields& fields() = 0;

    virtual void compute_source_term(const FieldState fstate) = 0;

    virtual void compute_nueff(const FieldState fstate) = 0;

    virtual void compute_diffusion_term(const FieldState fstate) = 0;

    virtual void compute_advection_term(const FieldState fstate) = 0;

    virtual void solve(const FieldState fstate) = 0;

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
    {}

    PDEFields& fields() override { return m_fields; }

    void compute_source_term(const FieldState fstate) override
    {
        m_src_op(fstate);
    }

    void compute_nueff(const FieldState ) override
    {
        // TBD
    }

    void compute_diffusion_term(const FieldState ) override
    {
        // m_diff_op(fstate);
    }

    void compute_advection_term(const FieldState ) override
    {
        // m_adv_op(fstate);
    }

    void solve(const FieldState ) override
    {
        // m_linsys_solver(fstate);
    }

protected:
    const SimTime& m_time;

    FieldRepo& m_repo;

    PDEFields m_fields;

    SrcTermOp<PDE> m_src_op;
    AdvectionOp<PDE, Scheme> m_adv_op;
    DiffusionOp<PDE> m_diff_op;
    LinSysSolve<PDE> m_linsys_solver;
};

} // namespace pde
} // namespace amr_wind

#endif /* PDE_H */
