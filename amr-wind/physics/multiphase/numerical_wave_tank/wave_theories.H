#ifndef WAVES_THEORIES_H_
#define WAVES_THEORIES_H_

#include <AMReX_FArrayBox.H>
#include <cmath>

namespace amr_wind {
namespace nwt {

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Vector<amrex::Real>
linear_monochromatic_waves(
    amrex::Real wavelength,
    amrex::Real waterdepth,
    amrex::Real amplitude,
    amrex::Real x,
    amrex::Real y,
    amrex::Real z,
    amrex::Real time)
{

    const amrex::Real k_ = 2.0 * M_PI / wavelength;
    const amrex::Real omega =
        std::pow(k_ * 9.81 * std::tanh(k_ * waterdepth), 0.5);
    const amrex::Real phase = k_ * x + 0. * y - omega * time;

    const amrex::Real eta = amplitude * std::sin(phase);
    const amrex::Real u_w = omega * amplitude *
                            std::cosh(k_ * (z + waterdepth)) /
                            std::sinh(k_ * waterdepth) * std::cos(phase);
    const amrex::Real v_w = 0.0;
    const amrex::Real w_w = omega * amplitude *
                            std::sinh(k_ * (z + waterdepth)) /
                            std::sinh(k_ * waterdepth) * std::sin(phase);
    const amrex::Vector<amrex::Real> out = {{eta, u_w, v_w, w_w}};
    return out;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Vector<amrex::Real>
stokes_waves(
    amrex::Real wavelength,
    amrex::Real waterdepth,
    amrex::Real amplitude,
    amrex::Real x,
    amrex::Real y,
    amrex::Real z,
    amrex::Real time,
    int order)
{
    const amrex::Real k_ = 2.0 * M_PI / wavelength;
    const amrex::Real omega =
        std::pow(k_ * 9.81 * std::tanh(k_ * waterdepth), 0.5);
    const amrex::Real phase = k_ * x + 0. * y - omega * time;

    const amrex::Real eta = 0.0;
    const amrex::Real u_w = 0.0;
    const amrex::Real v_w = 0.0;
    const amrex::Real w_w = 0.0;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
stokes_coefficients(amrex::Real k, amrex::Real waterdepth, int StokesOrder)
{
    amrex::Real kd = k * waterdepth;
    if (kd > 50 * M_PI) kd = 50. * M_PI; // Limited value
}

} // namespace nwt
} // namespace amr_wind

#endif