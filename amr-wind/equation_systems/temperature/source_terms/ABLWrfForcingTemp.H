#ifndef ABLWRFFORCINGTEMP_H
#define ABLWRFFORCINGTEMP_H

#include <AMReX_REAL.H>
#include <memory>
#include "amr-wind/core/FieldRepo.H"
#include "amr-wind/equation_systems/temperature/TemperatureSource.H"
#include "amr-wind/core/SimTime.H"
#include "amr-wind/utilities/FieldPlaneAveraging.H"
#include "amr-wind/utilities/ncutils/nc_interface.H"
#include "amr-wind/wind_energy/ABLMesoscaleForcing.H"
#include "amr-wind/wind_energy/ABLMesoscaleInput.H"

namespace amr_wind {
namespace pde {
namespace temperature {

class ABLWrfForcingTemp
    : public ABLMesoscaleForcing
    , public TemperatureSource::Register<ABLWrfForcingTemp>
{
public:
    static std::string identifier() { return "ABLWrfForcingTemp"; }

    explicit ABLWrfForcingTemp(const CFDSim& sim);

    ~ABLWrfForcingTemp() override;

    void operator()(
        const int lev,
        const amrex::MFIter& mfi,
        const amrex::Box& bx,
        const FieldState fstate,
        const amrex::Array4<amrex::Real>& src_term) const override;

    void mean_temperature_init(const ABLMesoscaleInput& ncfile);

    void mean_temperature_init(
        const FieldPlaneAveraging& tavg, const ABLMesoscaleInput& ncfile);

    amrex::Real mean_temperature_heights(std::unique_ptr<ABLMesoscaleInput>& ncfile);

    amrex::Real mean_temperature_heights(
        const FieldPlaneAveraging& tavg, std::unique_ptr<ABLMesoscaleInput>& ncfile);

    amrex::Vector<amrex::Real>& theta_wrf_error() { return m_err_Theta; }

private:
    // these are inputs
    amrex::Gpu::DeviceVector<amrex::Real> m_wrf_ht;
    amrex::Gpu::DeviceVector<amrex::Real> m_wrf_theta_vals;

    // this is the instantaneous planar average (at AMR-Wind levels)
    amrex::Gpu::DeviceVector<amrex::Real> m_theta_ht;
    amrex::Gpu::DeviceVector<amrex::Real> m_theta_vals;

    // this specifies the source term
    amrex::Gpu::DeviceVector<amrex::Real> m_error_wrf_avg_theta;

    amrex::Vector<amrex::Real> m_err_Theta;

    // this should be m_ind_polyOrder+1 (currently, cubic polynomial is assumed)
    amrex::Array<amrex::Real, 4> m_poly_coeff_theta;
};

} // namespace temperature

} // namespace pde

} // namespace amr_wind

#endif
