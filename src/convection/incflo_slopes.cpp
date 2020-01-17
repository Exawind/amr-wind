#include <incflo.H>

using namespace amrex;

//
// Compute the slopes of Sborder (velocity, density or tracer)
//
void
incflo::incflo_compute_slopes (int lev, Real time, MultiFab& Sborder,
                               Vector<std::unique_ptr<MultiFab>>& xslopes_in,
                               Vector<std::unique_ptr<MultiFab>>& yslopes_in,
                               Vector<std::unique_ptr<MultiFab>>& zslopes_in,
                               int slopes_comp, int ncomp)
{
}
