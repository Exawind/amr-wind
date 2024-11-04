#ifndef TEMPSPONGEFORCING_H
#define TEMPSPONGEFORCING_H

#include "amr-wind/equation_systems/temperature/TemperatureSource.H"
#include "amr-wind/core/SimTime.H"
#include "amr-wind/CFDSim.H"

namespace amr_wind::pde::temperature {

class TempSpongeForcing : public TemperatureSource::Register<TempSpongeForcing>
{
public:
    static std::string identifier() { return "TempSpongeForcing"; }

    explicit TempSpongeForcing(const CFDSim& sim);

    ~TempSpongeForcing() override;

    void operator()(
        const int lev,
        const amrex::MFIter& mfi,
        const amrex::Box& bx,
        const FieldState /*fstate*/,
        const amrex::Array4<amrex::Real>& src_term) const override;

private:
    const amrex::AmrCore& m_mesh;
    const Field& m_temperature;
    amrex::Vector<amrex::Real> m_theta_heights;
    amrex::Vector<amrex::Real> m_theta_values;
    amrex::Gpu::DeviceVector<amrex::Real> m_theta_heights_d;
    amrex::Gpu::DeviceVector<amrex::Real> m_theta_values_d;
    amrex::Real m_sponge_start{600};
};

} // namespace amr_wind::pde::temperature
#endif