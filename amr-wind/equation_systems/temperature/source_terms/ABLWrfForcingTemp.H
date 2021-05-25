#ifndef ABLWRFFORCINGTEMP_H
#define ABLWRFFORCINGTEMP_H

#include <AMReX_REAL.H>
#include <AMReX_Vector.H>
#include <memory>
#include <sys/_types/_size_t.h>
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

    amrex::Real mean_temperature_heights(
        const FieldPlaneAveraging&, std::unique_ptr<ABLWRFfile>&);

    amrex::Real mean_temperature_heights(
        const amrex::Vector<amrex::Real>&, std::unique_ptr<ABLWRFfile>&);

    void mean_temperature_init(const FieldPlaneAveraging&, const ABLWRFfile&);

    void invertMat(const amrex::Array2D<amrex::Real,0,3,0,3>&, amrex::Array2D<amrex::Real,0,3,0,3>&);

    void indirectForcingInit();

private:
    const SimTime& m_time;
    const amrex::AmrCore& m_mesh;
    std::string m_forcing_scheme;

    amrex::Gpu::DeviceVector<amrex::Real> m_theta_ht;
    amrex::Gpu::DeviceVector<amrex::Real> m_theta_vals;
    amrex::Gpu::DeviceVector<amrex::Real> m_wrf_theta_vals;
    amrex::Gpu::DeviceVector<amrex::Real> m_wrf_ht;

    amrex::Gpu::DeviceVector<amrex::Real> m_error_wrf_avg_theta; 

    int m_idx_time;

    int m_axis{2};

    amrex::Vector<amrex::Real> m_zht;

    int m_nht;
    size_t m_nlevels;

    amrex::Array2D<amrex::Real,0,3,0,3> m_im_zTz;

    amrex::Array<amrex::Real, 4> m_poly_coeff_theta;

    amrex::Real m_gain_coeff{0.2};

};

} // namespace temperature

} // namespace pde

} // namespace amr_wind

#endif
