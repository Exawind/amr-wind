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

    static void invertMat(
        const amrex::Array2D<amrex::Real, 0, 3, 0, 3>& m,
        amrex::Array2D<amrex::Real, 0, 3, 0, 3>& im);

protected:
    const SimTime& m_time;
    const amrex::AmrCore& m_mesh;

    int m_idx_time;
    int m_axis{2};

    // inputs
    std::string m_forcing_scheme;
    amrex::Real m_gain_coeff{0.2}; // relaxation factor for all forcing schemes
    bool m_debug{false};           // print out profiles for debugging

    // planar averaging heights
    amrex::Vector<amrex::Real> m_zht;
    int m_nht;

    //
    // regression vars for indirect forcing
    //
    // int m_ind_polyOrder;
    int m_norm_zmax = 0;
    amrex::Real m_scaleFact =
        1e-3; // to improve conditioning of regression matrix

    amrex::Array2D<amrex::Real, 0, 3, 0, 3> m_im_zTz;

    // - input (or current value of, if forcing_transition is specified)
    // regression weighting profile
    std::string m_forcing_transition;
    amrex::Real m_transition_height; // constant or variable
    bool m_update_transition_height = false;
    amrex::Real m_transition_thickness; // constant

    // - WORKAROUND
    amrex::Vector<amrex::Real> m_transition_height_hist;

    amrex::Vector<amrex::Real> m_weighting_heights; // constant or variable
    amrex::Vector<amrex::Real> m_weighting_values;  // constant or variable

    // - W(z) at planar averaging heights
    amrex::Vector<amrex::Real> m_W;

    //
    // weighting profile manipulation
    //
    void setTransitionWeighting(); // will override weighting profile
    void
    updateWeights(); // interpolate current profile to planar averaging heights

    //
    // forcing profile manipulation
    //
    bool forcingToConstant()
    {
        // - check if m_forcing_transition ends with "Constant"
        if (m_forcing_transition.length() < 8) {
            return false;
        }
        return (
            0 == m_forcing_transition.compare(
                 m_forcing_transition.length() - 8, 8, "Constant"));
    }

    // - slope varies linearly from the actual gradient of the forcing profile
    // at m_transition_height
    //   to 0 at (m_transition_height + m_transition_thicnkess)
    void constantForcingTransition(amrex::Vector<amrex::Real>& error);

    // - blend between a lower and upper forcing profile; the blending fraction
    // is assumed to be
    //   equal to W(z), but this is not a requirement; W=1 and W=0 correspond to
    //   the lower and upper profiles, respectively
    void blendForcings(
        const amrex::Vector<amrex::Real> lower,
        const amrex::Vector<amrex::Real> upper,
        amrex::Vector<amrex::Real>& error);
};

class ABLWRFfile
{

public:
    ABLWRFfile(const std::string wrfFile);

    // ~ABLWRFfile() = default;

    const amrex::Vector<amrex::Real>& wrf_heights() const;

    const amrex::Vector<amrex::Real>& wrf_times() const;

    const amrex::Vector<amrex::Real>& wrf_u() const;

    const amrex::Vector<amrex::Real>& wrf_v() const;

    const amrex::Vector<amrex::Real>& wrf_temp() const;

    const amrex::Vector<amrex::Real>& wrf_tflux() const;

    // ***FIXME***
    // const amrex::Vector<amrex::Real>& wrf_transition_height() const;

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
