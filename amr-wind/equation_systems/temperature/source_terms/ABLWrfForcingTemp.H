#ifndef ABLWRFFORCINGTEMP_H
#define ABLWRFFORCINGTEMP_H

#include "amr-wind/equation_systems/temperature/TemperatureSource.H"

namespace amr_wind {
namespace pde {
namespace temperature {

class ABLWrfForcingTemp : public TemperatureSource::Register<ABLWrfForcingTemp>
{
public:
    static const std::string identifier() { return "ABLWrfForcingTemp"; }

    explicit ABLWrfForcingTemp(const CFDSim& sim);
    virtual ~ABLWrfForcingTemp();
    virtual void operator()(
        const int lev,
        const amrex::MFIter& mfi,
        const amrex::Box& bx,
        const FieldState fstate,
        const amrex::Array4<amrex::Real>& src_term) const override;

    void read_forcing_file();

    void mean_temperature_heights(const FieldPlaneAveraging&);

    void mean_temperature_init(const FieldPlaneAveraging&);
};

} // namespace temperature

} // namespace pde

} // namespace amr_wind

#endif
