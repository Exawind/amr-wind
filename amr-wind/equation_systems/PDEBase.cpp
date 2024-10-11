#include "amr-wind/equation_systems/PDEBase.H"
#include "amr-wind/equation_systems/SchemeTraits.H"
#include "amr-wind/CFDSim.H"
#include "amr-wind/core/FieldRepo.H"
#include "amr-wind/equation_systems/PDEHelpers.H"
#include "amr-wind/incflo_enums.H"
#include "amr-wind/core/field_ops.H"

#include "AMReX_ParmParse.H"

namespace amr_wind::pde {

PDEFields::PDEFields(FieldRepo& repo_in, const std::string& var_name)
    : repo(repo_in)
    , field(repo.get_field(var_name))
    , mueff(repo.get_field(pde_impl::mueff_name(var_name)))
    , src_term(repo.get_field(pde_impl::src_term_name(var_name)))
    , diff_term(repo.get_field(pde_impl::diff_term_name(var_name)))
    , conv_term(repo.get_field(pde_impl::conv_term_name(var_name)))
{}

PDEMgr::PDEMgr(CFDSim& sim) : m_sim(sim)
{
    amrex::ParmParse pp("incflo");
    pp.query("use_godunov", m_use_godunov);
    pp.query("constant_density", m_constant_density);

    bool is_anelastic = false;
    {
        amrex::ParmParse pp_abl("ABL");
        pp_abl.query("anelastic", is_anelastic);
    }

    if ((is_anelastic) && (m_constant_density)) {
        amrex::Print()
            << "WARNING: Anelastic implies variable density. Setting "
               "contant_density to false. Add incflo.constant_density=false to "
               "your input file to remove this warning."
            << std::endl;
        m_constant_density = false;
    }

    m_scheme =
        m_use_godunov ? fvm::Godunov::scheme_name() : fvm::MOL::scheme_name();
}

PDEBase& PDEMgr::register_icns()
{
    const std::string name = "ICNS-" + m_scheme;

    m_icns = amr_wind::pde::PDEBase::create(name, m_sim);
    return *m_icns;
}

PDEBase& PDEMgr::register_transport_pde(const std::string& pde_name)
{
    const std::string name = pde_name + "-" + m_scheme;

    if (contains(name)) {
        amrex::Print() << "WARNING: Multiple requests to register PDE: "
                       << pde_name << std::endl;
        return operator()(name);
    }

    return create(name, m_sim);
}

bool PDEMgr::has_pde(const std::string& pde_name) const
{
    const std::string name = pde_name + "-" + m_scheme;
    return contains(name);
}

int PDEMgr::num_ghost_state() const
{
    return m_use_godunov ? fvm::Godunov::nghost_state : fvm::MOL::nghost_state;
}

void PDEMgr::advance_states()
{
    if (m_constant_density) {
        m_sim.repo().get_field("density").advance_states();
    }

    icns().fields().field.advance_states();
    for (auto& eqn : scalar_eqns()) {
        eqn->fields().field.advance_states();
    }
}

void PDEMgr::prepare_boundaries()
{
    // If state variables exist at NPH, fill their boundary cells
    const auto nph_time =
        m_sim.time().current_time() + 0.5 * m_sim.time().delta_t();
    if (m_constant_density &&
        m_sim.repo().field_exists("density", FieldState::NPH)) {
        auto& nph_field =
            m_sim.repo().get_field("density").state(FieldState::NPH);
        auto& new_field =
            m_sim.repo().get_field("density").state(FieldState::New);
        // Need to check if this step is necessary
        // Guessing that for no-op dirichlet bcs, the new needs to be copied to
        // NPH
        field_ops::copy(
            nph_field, new_field, 0, 0, new_field.num_comp(),
            new_field.num_grow());
        // Need to check that this uses the correct fop and fills cells at this
        // state
        nph_field.fillphysbc(nph_time);
    }

    if (m_sim.repo().field_exists(
            icns().fields().field.base_name(), FieldState::NPH)) {
        auto& nph_field = icns().fields().field.state(FieldState::NPH);
        auto& new_field = icns().fields().field.state(FieldState::New);
        field_ops::copy(
            nph_field, new_field, 0, 0, new_field.num_comp(),
            new_field.num_grow());
        nph_field.fillphysbc(nph_time);
    }
    for (auto& eqn : scalar_eqns()) {
        eqn->fields().field.advance_states();
        if (m_sim.repo().field_exists(
                eqn->fields().field.base_name(), FieldState::NPH)) {
            auto& nph_field = eqn->fields().field.state(FieldState::NPH);
            auto& new_field = eqn->fields().field.state(FieldState::New);
            field_ops::copy(
                nph_field, new_field, 0, 0, new_field.num_comp(),
                new_field.num_grow());
            nph_field.fillphysbc(nph_time);
        }
    }
}

void PDEMgr::fillpatch_state_fields(
    const amrex::Real time, const FieldState fstate)
{
    if (m_constant_density) {
        m_sim.repo().get_field("density").state(fstate).fillpatch(time);
    }

    icns().fields().field.state(fstate).fillpatch(time);
    for (auto& eqn : scalar_eqns()) {
        eqn->fields().field.state(fstate).fillpatch(time);
    }
}

} // namespace amr_wind::pde
