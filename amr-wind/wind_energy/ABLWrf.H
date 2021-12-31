#ifndef ABLWRF_H
#define ABLWRF_H

#include <string>
#include "AMReX_REAL.H"
#include "AMReX_Array.H"
#include "amr-wind/CFDSim.H"

namespace amr_wind {

class ABLWrfForcing
{

public:
    ABLWrfForcing(const CFDSim& sim, const std::string identifier);

    void indirectForcingInit();
    void invertMat(
        const amrex::Array2D<amrex::Real,0,3,0,3>&,
        amrex::Array2D<amrex::Real,0,3,0,3>&);

protected:
    const SimTime& m_time;
    const amrex::AmrCore& m_mesh;

    std::string m_forcing_scheme;
    amrex::Real m_gain_coeff{0.2};
    int m_norm_zmax = 0;
    amrex::Real m_scaleFact = 1e-3; // to improve conditioning of regression matrix
    int m_idx_time;
    int m_axis{2};

    // planar averaging heights
    amrex::Vector<amrex::Real> m_zht;
    int m_nht;

    // regression vars for indirect forcing
  //int m_ind_polyOrder;
    amrex::Array2D<amrex::Real,0,3,0,3> m_im_zTz;

    // - input (or current, if forcing_transition is specified) regression weighting profile
    std::string m_forcing_transition;
    amrex::Real m_transition_height; // constant or variable
    bool m_update_transition_height = false;
    amrex::Real m_transition_thickness; // constant
    amrex::Vector<amrex::Real> m_weighting_heights; // constant or variable
    amrex::Vector<amrex::Real> m_weighting_values; // constant or variable

    // - W(z) at planar averaging heights
    amrex::Vector<amrex::Real> m_W;

    // interpolate current weighting profile to planar averaging heights
    void setTransitionWeighting(); // will override weighting profile
    void updateWeights(); // interpolate current profile to planar averaging heights

    //
    // forcing profile manipulation
    //
    bool forcingToConstant() {
        // - check if m_forcing_transition ends with "Constant"
        if (m_forcing_transition.length() >= 8) {
            return (0 == m_forcing_transition.compare(m_forcing_transition.length()-8, 8, "Constant"));
        } else {
            return false;
        }
    }

    // - slope varies linearly from the actual gradient of the forcing profile at m_transition_height
    //   to 0 at (m_transition_height + m_transition_thicnkess)
    void constantForcingTransition(amrex::Vector<amrex::Real>&);

};

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

    bool is_wrf_tendency_forcing() const;

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
    bool m_abl_wrf_tendency;
};

} // namespace amr_wind

#endif
