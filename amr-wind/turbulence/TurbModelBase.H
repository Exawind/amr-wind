#ifndef TURBMODELBASE_H
#define TURBMODELBASE_H

#include "amr-wind/turbulence/TurbModel.H"
#include "amr-wind/utilities/IOManager.H"
#include "amr-wind/core/field_ops.H"

namespace amr_wind::turbulence {

namespace turb_base_impl {

// For transport model with constant properties implement specializations that
// avoid creation of an intermediate scratch buffer.

template <
    typename Transport,
    typename std::enable_if<Transport::constant_properties>::type* = nullptr>
inline void visc_update(Field& evisc, Field& tvisc, Transport& transport)
{
    evisc.setVal(transport.viscosity());
    field_ops::saxpy(
        evisc, 1.0, tvisc, 0, 0, evisc.num_comp(), evisc.num_grow());
}

template <
    typename Transport,
    typename std::enable_if<Transport::constant_properties>::type* = nullptr>
inline void alpha_update(Field& evisc, Field& tvisc, Transport& transport)
{
    evisc.setVal(transport.thermal_diffusivity());
    field_ops::saxpy(
        evisc, 1.0 / transport.turbulent_prandtl(), tvisc, 0, 0,
        evisc.num_comp(), evisc.num_grow());
}

template <
    typename Transport,
    typename std::enable_if<Transport::constant_properties>::type* = nullptr>
inline void scal_diff_update(
    Field& evisc, Field& tvisc, Transport& transport, const std::string& name)
{
    evisc.setVal(transport.viscosity() / transport.laminar_schmidt(name));
    field_ops::saxpy(
        evisc, 1.0 / transport.turbulent_schmidt(name), tvisc, 0, 0,
        evisc.num_comp(), evisc.num_grow());
}

template <
    typename Transport,
    typename std::enable_if<!Transport::constant_properties>::type* = nullptr>
inline void visc_update(Field& evisc, Field& tvisc, Transport& transport)
{
    auto lam_mu = transport.mu();
    field_ops::lincomb(
        evisc, 1.0, *lam_mu, 0, 1.0, tvisc, 0, 0, evisc.num_comp(),
        evisc.num_grow());
}

template <
    typename Transport,
    typename std::enable_if<!Transport::constant_properties>::type* = nullptr>
inline void alpha_update(Field& evisc, Field& tvisc, Transport& transport)
{
    auto lam_alpha = transport.alpha();
    field_ops::lincomb(
        evisc, 1.0, *lam_alpha, 0, 1.0 / transport.turbulent_prandtl(), tvisc,
        0, 0, evisc.num_comp(), evisc.num_grow());
}

template <
    typename Transport,
    typename std::enable_if<!Transport::constant_properties>::type* = nullptr>
inline void scal_diff_update(
    Field& evisc, Field& tvisc, Transport& transport, const std::string& name)
{
    auto lam_mu = transport.mu();
    field_ops::lincomb(
        evisc, 1.0 / transport.laminar_schmidt(name), *lam_mu, 0,
        1.0 / transport.turbulent_schmidt(name), tvisc, 0, 0, evisc.num_comp(),
        evisc.num_grow());
}

} // namespace turb_base_impl

/** Turbulence model aspects common to both LES and RANS models
 *  \ingroup turb_model
 *
 *  This class implements more of the common methods from TurbulenceModel that
 *  are common to most turbulence model. This is different from TurbModel in
 *  that that class is also valid for Laminar which isn't really a turbulence
 *  model but is a convenient way to represent the absence of a turbulence
 *  model.
 */
template <typename Transport>
class TurbModelBase : public TurbModel<Transport>
{
public:
    explicit TurbModelBase(CFDSim& sim)
        : TurbModel<Transport>(sim)
        , m_mu_turb(sim.repo().declare_field("mu_turb", 1, 1, 1))
    {
        this->m_sim.io_manager().register_io_var("mu_turb");

        // Effective viscosity field
        auto& fields = this->m_sim.pde_manager().icns().fields();
        this->m_sim.io_manager().register_io_var(fields.mueff.name());

        m_mu_turb.set_default_fillpatch_bc(this->m_sim.time());
    }

    //! Return the turbulent viscosity field
    Field& mu_turb() override { return m_mu_turb; }

    //! Return the turbulent thermal diffusivity field
    Field& alpha_turb() override
    {
        AMREX_ASSERT(m_alpha_turb != nullptr);
        return *m_alpha_turb;
    }

    // clang-format off

    /** Update the effective viscosity field
     *
     *  \f$\mu_\mathrm{eff} = \mu + \mu_t\f$
     *
     *  \param mueff Effective viscosity field
     */
    void update_mueff(Field& mueff) override
    {
        turb_base_impl::visc_update(mueff, this->m_mu_turb, this->m_transport);
    }

    /** Update the effective thermal diffusivity field
     *
     *  \f$\alpha_\mathrm{eff} = \alpha + \alpha_t\f$
     *
     *  \param mueff Effective viscosity field
     */
    void update_alphaeff(Field& alphaeff) override
    {
        turb_base_impl::alpha_update(alphaeff, this->m_mu_turb, this->m_transport);
    }

    //! Interface to update scalar diffusivity based on Schmidt number
    void update_scalar_diff(Field& deff, const std::string& name) override
    {
        turb_base_impl::scal_diff_update(deff, this->m_mu_turb, this->m_transport, name);
    }

    // clang-format on

protected:
    //! Reference to the turbulent viscosity field
    Field& m_mu_turb;

    //! Reference to thermal diffusivity field
    Field* m_alpha_turb{nullptr};
};

} // namespace amr_wind::turbulence

#endif /* TURBMODELBASE_H */
