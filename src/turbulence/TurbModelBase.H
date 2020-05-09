#ifndef TURBMODELBASE_H
#define TURBMODELBASE_H

#include "TurbModel.H"
#include "IOManager.H"
#include "field_ops.H"

namespace amr_wind {
namespace turbulence {

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
        evisc, 1.0, lam_alpha, 0, 1.0 / transport.turbulent_prandtl(),
        tvisc, 0, 0, evisc.num_comp(), evisc.num_grow());
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
        1.0 / transport.turbulent_schmidt(), tvisc, 0, 0, evisc.num_comp(),
        evisc.num_grow());
}

} // namespace turb_base_impl

/** Turbulence model aspects common to both LES and RANS models
 */
template<typename Transport>
class TurbModelBase : public TurbModel<Transport>
{
public:
    TurbModelBase(CFDSim& sim)
        : TurbModel<Transport>(sim)
        , m_mu_turb(sim.repo().declare_field("mu_turb", 1, 1, 1))
    {
        this->m_sim.io_manager().register_io_var("mu_turb");

        // Effective viscosity field
        auto& fields = this->m_sim.pde_manager().icns().fields();
        this->m_sim.io_manager().register_io_var(fields.mueff.name());
    }

    virtual Field& mu_turb() override { return m_mu_turb; }

    virtual Field& alpha_turb() override
    {
        AMREX_ASSERT(m_alpha_turb != nullptr);
        return *m_alpha_turb;
    }

    // clang-format off

    //! Interface to update effective viscosity (mu_eff = mu + mu_t)
    virtual void update_mueff(Field& mueff) override
    {
        turb_base_impl::visc_update(mueff, this->m_mu_turb, this->m_transport);
    }

    //! Interface to update effective thermal diffusivity
    virtual void update_alphaeff(Field& alphaeff) override
    {
        turb_base_impl::alpha_update(alphaeff, this->m_mu_turb, this->m_transport);
    }

    //! Interface to update scalar diffusivity based on Schmidt number
    virtual void update_scalar_diff(Field& deff, const std::string& name) override
    {
        turb_base_impl::scal_diff_update(deff, this->m_mu_turb, this->m_transport, name);
    }

    // clang-format on

protected:
    Field& m_mu_turb;

    Field* m_alpha_turb{nullptr};
};

}
}

#endif /* TURBMODELBASE_H */
