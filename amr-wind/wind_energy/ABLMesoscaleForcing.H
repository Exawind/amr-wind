#ifndef ABLMESOSCALEFORCING_H
#define ABLMESOSCALEFORCING_H

#include <string>
#include "AMReX_REAL.H"
#include "AMReX_Array.H"
#include "AMReX_Vector.H"
#include "amr-wind/CFDSim.H"

namespace amr_wind {

class ABLMesoscaleForcing
{

public:
    ABLMesoscaleForcing(const CFDSim& sim, const std::string& identifier);

    void indirectForcingInit();

    static void invertMat(
        const amrex::Array2D<amrex::Real, 0, 3, 0, 3>& m,
        amrex::Array2D<amrex::Real, 0, 3, 0, 3>& im);

protected:
    const SimTime& m_time;
    const amrex::AmrCore& m_mesh;

    int m_idx_time{0};
    int m_axis{2};

    // inputs
    std::string m_forcing_scheme;
    amrex::Real m_gain_coeff{0.2}; // relaxation factor for all forcing schemes
    bool m_debug{false};           // print out profiles for debugging

    // planar averaging heights
    amrex::Vector<amrex::Real> m_zht;
    int m_nht{0}; // set by derived ABLMesoForcing[Mom|Temp] class

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
    bool m_update_transition_height{false};
    amrex::Real m_transition_thickness; // constant

    // - WORKAROUND
    amrex::Vector<amrex::Real> m_transition_height_hist;

    amrex::Vector<amrex::Real> m_weighting_heights; // constant or variable
    amrex::Vector<amrex::Real> m_weighting_values;  // constant or variable
    bool m_user_specified_weighting{false};

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

} // namespace amr_wind

#endif
