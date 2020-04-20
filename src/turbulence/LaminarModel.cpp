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
laminar_visc_update(Field& evisc, ScratchField&, const Transport& transport)
{
    evisc.setVal(transport.viscosity());
}

template <
    typename Transport,
    typename std::enable_if<!Transport::constant_properties>::type* = nullptr>
inline void laminar_visc_update(Field& evisc, ScratchField& visc)
{
    field_ops::copy(evisc, visc, 0, 0, evisc.num_comp(), evisc.num_grow());
}

template <
    typename Transport,
    typename std::enable_if<Transport::constant_properties>::type* = nullptr>
inline void
laminar_alpha_update(Field& evisc, ScratchField&, const Transport& transport)
{
    evisc.setVal(transport.thermal_diffusivity());
}

template <
    typename Transport,
    typename std::enable_if<!Transport::constant_properties>::type* = nullptr>
inline void laminar_alpha_update(Field& evisc, ScratchField& visc)
{
    field_ops::copy(evisc, visc, 0, 0, evisc.num_comp(), evisc.num_grow());
}

} // namespace

template <typename Transport>
void Laminar<Transport>::update_turbulent_viscosity()
{
    AMREX_ASSERT(this->m_mueff != nullptr);
    laminar_visc_update(*this->m_mueff, *this->mu(), this->m_transport);

    if (this->m_alphaeff != nullptr)
        laminar_alpha_update(
            *this->m_alphaeff, *this->alpha(), this->m_transport);
}

} // namespace turbulence

INSTANTIATE_TURBULENCE_MODEL(Laminar);

} // namespace amr_wind
