#ifndef WAVES_THEORIES_H_
#define WAVES_THEORIES_H_

#include <AMReX_FArrayBox.H>
#include <cmath>

namespace amr_wind {
namespace nwt {

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void linear_monochromatic_waves(
    amrex::Real wavelength,
    amrex::Real waterdepth,
    amrex::Real amplitude,
    amrex::Real x,
    amrex::Real y,
    amrex::Real z,
    amrex::Real time,
    amrex::Real& eta,
    amrex::Real& u_w,
    amrex::Real& v_w,
    amrex::Real& w_w)
{

    const amrex::Real k_ = 2.0 * M_PI / wavelength;
    const amrex::Real omega =
        std::pow(k_ * 9.81 * std::tanh(k_ * waterdepth), 0.5);
    const amrex::Real phase = k_ * x + 0. * y - omega * time;

    eta = amplitude * std::cos(phase);
    u_w = omega * amplitude * std::cosh(k_ * (z + waterdepth)) /
          std::sinh(k_ * waterdepth) * std::cos(phase);
    v_w = 0.0;
    w_w = omega * amplitude * std::sinh(k_ * (z + waterdepth)) /
          std::sinh(k_ * waterdepth) * std::sin(phase);
}

} // namespace nwt
} // namespace amr_wind

#endif
