#include <incflo.H>

#ifdef AMREX_USE_EB
#include <AMReX_EB_utils.H>
#endif

using namespace amrex;

//
// Compute the MAC (normal) velocities on face centroids
//
void
incflo::incflo_compute_MAC_velocities( Vector< std::unique_ptr<MultiFab> >& vel_forces_in,
                                       Vector< std::unique_ptr<MultiFab> >& scal_forces_in,
                                       Vector< std::unique_ptr<MultiFab> >& vel_in,
                                       Vector< std::unique_ptr<MultiFab> >& density_in,
                                       Vector< std::unique_ptr<MultiFab> >& tracer_in,
                                       Real time)
{
}
