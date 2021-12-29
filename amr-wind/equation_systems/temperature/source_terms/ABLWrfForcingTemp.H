#ifndef ABLWRFFORCINGTEMP_H
#define ABLWRFFORCINGTEMP_H

#include <AMReX_REAL.H>
#include <memory>
#include "amr-wind/core/FieldRepo.H"
#include "amr-wind/equation_systems/temperature/TemperatureSource.H"
#include "amr-wind/core/SimTime.H"
#include "amr-wind/utilities/FieldPlaneAveraging.H"
#include "amr-wind/utilities/ncutils/nc_interface.H"
#include "amr-wind/wind_energy/ABLWrf.H"

namespace amr_wind {
namespace pde {
namespace temperature {

class ABLWrfForcingTemp : public ABLWrfForcing,
                          public TemperatureSource::Register<ABLWrfForcingTemp>
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

    amrex::Real mean_temperature_heights(
        const FieldPlaneAveraging&, std::unique_ptr<ABLWRFfile>&);

    amrex::Real mean_temperature_heights(std::unique_ptr<ABLWRFfile>&);
  
    void mean_temperature_init(const FieldPlaneAveraging&, const ABLWRFfile&);

    void mean_temperature_init(const ABLWRFfile&);

  amrex::Vector<amrex::Real>& theta_wrf_error(){
    return m_err_Theta;} 

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
