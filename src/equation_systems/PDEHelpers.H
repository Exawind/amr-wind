#ifndef PDEHELPERS_H
#define PDEHELPERS_H

#include <string>

#include "SimTime.H"
#include "FieldRepo.H"
#include "FieldFillPatchOps.H"
#include "FieldBCOps.H"

namespace amr_wind {
namespace pde_impl {

//! Field name that holds the convective term
inline std::string conv_term_name(const std::string& var)
{
    return var + "_conv_term";
}

//! Field name that holds the diffusion term (divtau and laplacian for scalars)
inline std::string diff_term_name(const std::string& var)
{
    return var + "_diff_term";
}

//! Field name for all the source terms
inline std::string src_term_name(const std::string& var)
{
    return var + "_src_term";
}

//! Effective viscosity for the transport equation
inline std::string mueff_name(const std::string& var) { return var + "_mueff"; }

} // namespace pde_impl

namespace pde {

struct PDEFields
{
    PDEFields(FieldRepo& repo_in, const std::string& var_name)
        : repo(repo_in)
        , field(repo.get_field(var_name))
        , mueff(repo.get_field(pde_impl::mueff_name(var_name)))
        , src_term(repo.get_field(pde_impl::src_term_name(var_name)))
        , diff_term(repo.get_field(pde_impl::diff_term_name(var_name)))
        , conv_term(repo.get_field(pde_impl::conv_term_name(var_name)))
    {}

    FieldRepo& repo;
    Field& field;
    Field& mueff;

    Field& src_term;
    Field& diff_term;
    Field& conv_term;
};

template<typename PDE, typename Scheme>
PDEFields create_fields_instance(const SimTime& time, FieldRepo& repo, const int probtype)
{
    repo.declare_field(PDE::var_name(), PDE::ndim, Scheme::nghost_state, Scheme::num_states);
    repo.declare_field(pde_impl::mueff_name(PDE::var_name()), 1, 1, 1);
    repo.declare_field(
        pde_impl::src_term_name(PDE::var_name()), PDE::ndim,
        Scheme::nghost_src, 1);
    repo.declare_field(
        pde_impl::diff_term_name(PDE::var_name()), PDE::ndim, 0,
        Scheme::num_diff_states);
    repo.declare_field(
        pde_impl::conv_term_name(PDE::var_name()), PDE::ndim, 0,
        Scheme::num_conv_states);

    PDEFields fields(repo, PDE::var_name());
    fields.field.register_fill_patch_op<FieldFillPatchOps<FieldBCDirichlet>>(
        repo.mesh(), time, probtype);
    fields.src_term.register_fill_patch_op<FieldFillPatchOps<FieldBCNoOp>>(
        repo.mesh(), time, probtype, FieldInterpolator::PiecewiseConstant);

    fields.field.fillpatch_on_regrid() = true;

    if (PDE::need_nph_state)
        fields.field.create_state(FieldState::NPH);

    return fields;
}

} // namespace pde
} // namespace amr_wind


#endif /* PDEHELPERS_H */
