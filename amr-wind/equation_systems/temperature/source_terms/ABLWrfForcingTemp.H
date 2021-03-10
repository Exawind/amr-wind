#ifndef ABLWRFFORCINGTEMP_H
#define ABLWRFFORCINGTEMP_H

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

  void mean_temperature_heights(const FieldPlaneAveraging&, std::unique_ptr<ABLWRFfile>&);
  
  void mean_temperature_init(const FieldPlaneAveraging&, const ABLWRFfile&);

 private:
  const SimTime& m_time;
  const amrex::AmrCore& m_mesh;
  std::string m_forcing_scheme{"direct"};

  amrex::Gpu::DeviceVector<amrex::Real> m_theta_ht;
  amrex::Gpu::DeviceVector<amrex::Real> m_theta_vals;
  amrex::Gpu::DeviceVector<amrex::Real> m_wrf_theta_vals;
  amrex::Gpu::DeviceVector<amrex::Real> m_wrf_ht;
  
  int m_idx_time;

  int m_axis{2};

};

} // namespace temperature

} // namespace pde

} // namespace amr_wind

#endif
