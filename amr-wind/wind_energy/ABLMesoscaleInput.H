#ifndef ABLMESOSCALEINPUT_H
#define ABLMESOSCALEINPUT_H

#include <string>
#include "AMReX_REAL.H"
#include "AMReX_Vector.H"

namespace amr_wind {

class ABLMesoscaleInput
{

public:
    ABLMesoscaleInput(const std::string ncfile);

    // ~ABLMesoscaleInput() = default;

    const amrex::Vector<amrex::Real>& wrf_heights() const;

    const amrex::Vector<amrex::Real>& wrf_times() const;

    const amrex::Vector<amrex::Real>& wrf_u() const;

    const amrex::Vector<amrex::Real>& wrf_v() const;

    const amrex::Vector<amrex::Real>& wrf_temp() const;

    const amrex::Vector<amrex::Real>& wrf_tflux() const;

    // ***FIXME***
    // const amrex::Vector<amrex::Real>& wrf_transition_height() const;

    bool is_wrf_tendency_forcing() const { return m_abl_wrf_tendency; }

    int nheights() const { return m_wrf_nheight; }
    int times() const { return m_wrf_ntime; }

private:
    std::string m_wrf_filename;

    amrex::Vector<amrex::Real> m_wrf_height;
    amrex::Vector<amrex::Real> m_wrf_time;

    amrex::Vector<amrex::Real> m_wrf_u;
    amrex::Vector<amrex::Real> m_wrf_v;
    amrex::Vector<amrex::Real> m_wrf_temp;
    amrex::Vector<amrex::Real> m_wrf_tflux;

    // optional input for partial profile assimilation
    // ***FIXME***
    // unexpected behaviors, as described in
    // ec5eb95c6ca853ce0fea8488e3f2515a2d6374e7 MUST COMMENT THIS LINE OUT
    // (variable def) to consistently fix problem:
    //
    // amrex::Vector<amrex::Real> m_wrf_transition_height;

    int m_wrf_nheight;
    int m_wrf_ntime;
    bool m_abl_wrf_tendency;
};

} // namespace amr_wind

#endif
