#ifndef ABLMESOSCALEINPUT_H
#define ABLMESOSCALEINPUT_H

#include <string>
#include "AMReX_REAL.H"
#include "AMReX_Vector.H"

namespace amr_wind {

class ABLMesoscaleInput
{

public:
    // cppcheck-suppress noExplicitConstructor
    ABLMesoscaleInput(std::string ncfile, std::string var_prefix = "");

    ~ABLMesoscaleInput() = default;

    const amrex::Vector<amrex::Real>& meso_heights() const { return m_height; }

    const amrex::Vector<amrex::Real>& meso_times() const { return m_time; }

    const amrex::Vector<amrex::Real>& meso_u() const { return m_u; }

    const amrex::Vector<amrex::Real>& meso_v() const { return m_v; }

    const amrex::Vector<amrex::Real>& meso_temp() const { return m_temp; }

    const amrex::Vector<amrex::Real>& meso_tflux() const { return m_tflux; }

    const amrex::Vector<amrex::Real>& meso_transition_height() const
    {
        return m_transition_height;
    }

    bool is_tendency_forcing() const { return m_abl_tendency; }

    int nheights() const { return m_nheight; }
    int times() const { return m_ntime; }

private:
    std::string m_filename;
    std::string m_var_prefix; // to support legacy WRF meso input files

    amrex::Vector<amrex::Real> m_height;
    amrex::Vector<amrex::Real> m_time;

    amrex::Vector<amrex::Real> m_u;
    amrex::Vector<amrex::Real> m_v;
    amrex::Vector<amrex::Real> m_temp;
    amrex::Vector<amrex::Real> m_tflux;

    // possible unexpected behaviors, as described in
    // ec5eb95c6ca853ce0fea8488e3f2515a2d6374e7
    //
    // See also:
    // - wind_energy/ABLMesoscaleInput.cpp
    // - equation_systems/icns/source_terms/ABLMesoForcingMom.cpp
    // - equation_systems/temperature/source_terms/ABLMesoForcingTemp.cpp
    //
    // workaround is to not try to read in and update this transition height
    // variable on the fly
    amrex::Vector<amrex::Real> m_transition_height;

    int m_nheight;
    int m_ntime;
    bool m_abl_tendency;
};

} // namespace amr_wind

#endif
