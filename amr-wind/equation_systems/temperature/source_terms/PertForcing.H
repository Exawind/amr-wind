#ifndef PERTFORCING_H
#define PERTFORCING_H

#include "amr-wind/equation_systems/temperature/TemperatureSource.H"
#include "amr-wind/core/SimTime.H"
#include "amr-wind/CFDSim.H"
#include "amr-wind/transport_models/TransportModel.H"

namespace amr_wind::pde::temperature {

class PertForcing : public TemperatureSource::Register<PertForcing>
{
public:
    static std::string identifier() { return "PertForcing"; }

    explicit PertForcing(const CFDSim& sim);

    ~PertForcing() override;

    void operator()(
        const int lev,
        const amrex::MFIter& mfi,
        const amrex::Box& bx,
        const FieldState /*fstate*/,
        const amrex::Array4<amrex::Real>& src_term) const override;

private:
    const CFDSim& m_sim;
    const SimTime& m_time;
    const amrex::AmrCore& m_mesh;
    amrex::Real m_xstart;
    amrex::Real m_xend;
    amrex::Real m_ystart;
    amrex::Real m_yend;
    amrex::Real m_zstart;
    amrex::Real m_zend;
    int m_time_index{100};
    int m_start_level{0};
    amrex::Real m_pert_amplitude{1.0};
};

} // namespace amr_wind::pde::temperature
#endif
