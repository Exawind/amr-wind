#ifndef TURBMODEL_H
#define TURBMODEL_H

#include "amr-wind/turbulence/TurbulenceModel.H"
#include "amr-wind/CFDSim.H"

namespace amr_wind {
namespace turbulence {

template<typename Transport>
class TurbModel : public TurbulenceModel
{
public:
    using TransportType = Transport;

    TurbModel(CFDSim& sim)
        : m_sim(sim)
        , m_transport(sim.repo())
    {}

    virtual void register_mueff_field(Field& mueff) override
    {
        m_mueff = &mueff;
    }

    virtual void register_alphaeff_field(Field& alphaeff) override
    {
        m_alphaeff = &alphaeff;
    }

    //! Return the dynamic viscosity field
    virtual std::unique_ptr<ScratchField> mu() override
    { return m_transport.mu(); }

    //! Return the thermal diffusivity field
    virtual std::unique_ptr<ScratchField> alpha() override
    { return m_transport.alpha(); }

    //! Return the scalar diffusivity field
    virtual std::unique_ptr<ScratchField> scalar_diffusivity(const std::string& name) override
    {
        return m_transport.scalar_diffusivity(name);
    }

    virtual Field& mueff() override
    {
        AMREX_ASSERT(m_mueff != nullptr);
        return *m_mueff;
    }

    virtual Field& alphaeff() override
    {
        AMREX_ASSERT(m_alphaeff != nullptr);
        return *m_alphaeff;
    }

protected:
    CFDSim& m_sim;

    //! Transport properties instance
    Transport m_transport;

    //! Reference to the effective viscosity field
    Field* m_mueff{nullptr};

    //! Reference to the thermal diffusivity field
    Field* m_alphaeff{nullptr};
};

}
}

#endif /* TURBMODEL_H */
