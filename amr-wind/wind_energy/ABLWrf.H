#ifndef ABLWRF_H
#define ABLWRF_H

#include <string>
#include "AMReX_REAL.H"
#include "AMReX_Array.H"
#include "amr-wind/CFDSim.H"

namespace amr_wind {

class ABLWRFfile
{

public:
    ABLWRFfile(const std::string);

    // ~ABLWRFfile() = default;

    const amrex::Vector<amrex::Real>& wrf_heights() const;

    const amrex::Vector<amrex::Real>& wrf_times() const;

    const amrex::Vector<amrex::Real>& wrf_u() const;

    const amrex::Vector<amrex::Real>& wrf_v() const;

    const amrex::Vector<amrex::Real>& wrf_temp() const;

    const amrex::Vector<amrex::Real>& wrf_tflux() const;

    int nheights() const;
    int times() const;

private:
    std::string m_wrf_filename;

    amrex::Vector<amrex::Real> m_wrf_height;
    amrex::Vector<amrex::Real> m_wrf_time;

    amrex::Vector<amrex::Real> m_wrf_u;
    amrex::Vector<amrex::Real> m_wrf_v;
    amrex::Vector<amrex::Real> m_wrf_temp;
    amrex::Vector<amrex::Real> m_wrf_tflux;

    int m_wrf_nheight;
    int m_wrf_ntime;
};

} // namespace amr_wind

#endif
