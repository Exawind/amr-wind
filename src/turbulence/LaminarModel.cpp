#include "LaminarModel.H"
#include "TurbModelDefs.H"
#include "field_ops.H"

namespace amr_wind {
namespace turbulence {

namespace {

// For laminar models we just copy the dynamic viscosity and thermal diffusivity
// into the effective/turbulent viscosity fields to satisfy the API. However, we
// use template specializations on the transport models to bypass an extra field
// creation/copy overhead when the properties are constant. The functions below
// perform a direct update using setVal for constant properties, but a copy from
// laminar field to effective field for non-constant transport properties.

template <
    typename Transport,
    typename std::enable_if<Transport::constant_properties>::type* = nullptr>
inline void
laminar_visc_update(Field& evisc, Laminar<Transport>&, const Transport& transport)
{
    evisc.setVal(transport.viscosity());
}

template <
    typename Transport,
    typename std::enable_if<!Transport::constant_properties>::type* = nullptr>
inline void laminar_visc_update(Field& evisc, Laminar<Transport>& lam, const Transport&)
{
    field_ops::copy(evisc, lam.mu(), 0, 0, evisc.num_comp(), evisc.num_grow());
}

template <
    typename Transport,
    typename std::enable_if<Transport::constant_properties>::type* = nullptr>
inline void
laminar_alpha_update(Field& evisc, Laminar<Transport>&, const Transport& transport)
{
    evisc.setVal(transport.thermal_diffusivity());
}

template <
    typename Transport,
    typename std::enable_if<!Transport::constant_properties>::type* = nullptr>
inline void laminar_alpha_update(Field& evisc, Laminar<Transport>& lam, const Transport&)
{
    field_ops::copy(evisc, lam.alpha(), 0, 0, evisc.num_comp(), evisc.num_grow());
}

template <
    typename Transport,
    typename std::enable_if<Transport::constant_properties>::type* = nullptr>
inline void
laminar_scal_diff_update(Field& evisc, Laminar<Transport>&, const Transport& transport,
                         const std::string& name)
{
    evisc.setVal(transport.viscosity() / transport.laminar_schmidt(name));
}

template <
    typename Transport,
    typename std::enable_if<!Transport::constant_properties>::type* = nullptr>
inline void laminar_scal_diff_update(Field& evisc, Laminar<Transport>& lam, const Transport&,
                                     const std::string& name)
{
    field_ops::copy(
        evisc, lam.scalar_diffusivity(name), 0, 0, evisc.num_comp(),
        evisc.num_grow());
}

} // namespace

template <typename Transport>
void Laminar<Transport>::update_turbulent_viscosity(const FieldState /* fstate */)
{
    // Empty function as there is no turbulent field
}

template<typename Transport>
void Laminar<Transport>::update_mueff(Field& mueff)
{
    laminar_visc_update(mueff, *this, this->m_transport);
}

template<typename Transport>
void Laminar<Transport>::update_alphaeff(Field& alphaeff)
{
    laminar_alpha_update(alphaeff, *this, this->m_transport);
}

template<typename Transport>
void Laminar<Transport>::update_scalar_diff(Field& deff, const std::string& name)
{
    laminar_scal_diff_update(deff, *this, this->m_transport, name);
}

} // namespace turbulence

INSTANTIATE_TURBULENCE_MODEL(Laminar);

} // namespace amr_wind
