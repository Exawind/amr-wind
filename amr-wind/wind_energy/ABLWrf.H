#ifdef ABLWRF_H
#define ABLWRF_H

#include <memory>

namespace amr_wind {

class ABLWRFfile
{

public:
    explicit ABLWRFfile(const std::string filewrf);

    // ~ABLWRFfile() = default;

    const amrex::Vector<amrex::Real>& wrf_heights() const
    {
        return m_wrf_height;
    }

    const amrex::Vector<amrex::Real>& wrf_times() const { return m_wrf_time; }

    const amrex::Vector<amrex::Real>& wrf_u() const { return m_wrf_u; }

    const amrex::Vector<amrex::Real>& wrf_v() const { return m_wrf_v; }

    const amrex::Vector<amrex::Real>& wrf_temp() const { return m_wrf_temp; }

    const amrex::Vector<amrex::Real>& wrf_tflux() const { return m_wrf_tflux; }

    int nheights() const { return m_nheight; }
    int times() const { return m_ntime; }

private:
    std::string m_wrf_filename;

    amrex::Vector<amrex::Real> m_wrf_height;
    amrex::Vector<amrex::Real> m_wrf_time;

    amrex::Vector<amrex::Real> m_wrf_u;
    amrex::Vector<amrex::Real> m_wrf_v;
    amrex::Vector<amrex::Real> m_wrf_temp;
    amrex::Vector<amrex::Real> m_wrf_tflux;

    int m_nheight;
    int m_ntime;
};

} // namespace amr_wind

#endif
