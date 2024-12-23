#ifndef DRAGTEMPFORCING_H
#define DRAGTEMPFORCING_H

#include "amr-wind/equation_systems/temperature/TemperatureSource.H"
#include "amr-wind/core/SimTime.H"
#include "amr-wind/CFDSim.H"
#include "amr-wind/transport_models/TransportModel.H"

namespace amr_wind::pde::temperature {

class DragTempForcing : public TemperatureSource::Register<DragTempForcing>
{
public:
    static std::string identifier() { return "DragTempForcing"; }

    explicit DragTempForcing(const CFDSim& sim);

    ~DragTempForcing() override;

    void operator()(
        const int lev,
        const amrex::MFIter& mfi,
        const amrex::Box& bx,
        const FieldState /*fstate*/,
        const amrex::Array4<amrex::Real>& src_term) const override;

private:
    const CFDSim& m_sim;
    const amrex::AmrCore& m_mesh;
    const Field& m_velocity;
    const Field& m_temperature;
    amrex::Real m_drag_coefficient{1.0};

    //! Transport model
    const transport::TransportModel& m_transport;

    //! Reference temperature
    std::unique_ptr<ScratchField> m_ref_theta;
};

} // namespace amr_wind::pde::temperature
#endif
