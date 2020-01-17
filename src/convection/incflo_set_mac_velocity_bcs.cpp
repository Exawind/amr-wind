
#include <incflo.H>

using namespace amrex;

//
// Set the BCs for face-centroid-based velocity components only
// 
void
incflo::set_MAC_velocity_bcs ( int lev,
                               Vector< std::unique_ptr<MultiFab> >& u_mac,
                               Vector< std::unique_ptr<MultiFab> >& v_mac,
                               Vector< std::unique_ptr<MultiFab> >& w_mac,
                               Real time)
{
}
