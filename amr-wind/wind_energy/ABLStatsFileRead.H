#ifndef ABLSTATSFILEREAD_H
#define ABLSTATSFILEREAD_H

#include <AMReX_INT.H>
#include <AMReX_Vector.H>
#include <string>
#include "AMReX_REAL.H"
#include "AMReX_Array.H"
#include "amr-wind/CFDSim.H"

namespace amr_wind {

class ABLReadStats
{

public:
    ABLReadStats(const std::string);

    const amrex::Vector<amrex::Real>& stats_ustar() const;
    const amrex::Vector<amrex::Real>& stats_time() const;
    const amrex::Vector<amrex::Real>& stats_theta() const;

    amrex::Real interpUstarTime(const amrex::Real);

  void interpThetaTime(const amrex::Real, amrex::Array<amrex::Real, 4>&);
    

private:
    std::string m_stat_filename;

    size_t m_stats_nt_steps;
    amrex::Vector<amrex::Real> m_stats_time;
    amrex::Vector<amrex::Real> m_stats_ustar;
 
    amrex::Vector<amrex::Real> m_stats_theta;
    amrex::Vector<amrex::Real> m_stats_theta_1D;

    amrex::Vector<amrex::Real> m_stats_v;
    amrex::Vector<amrex::Real> m_stats_v_1D;

    amrex::Vector<amrex::Real> m_stats_u;
    amrex::Vector<amrex::Real> m_stats_u_1D;

    amrex::Vector<amrex::Real> m_stats_vmag;
    amrex::Vector<amrex::Real> m_stats_vmag_1D;

    size_t m_stats_nlevels;
};

} // namespace amr_wind

#endif
