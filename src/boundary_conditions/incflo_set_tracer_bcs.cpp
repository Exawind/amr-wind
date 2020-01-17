#include <incflo.H>

using namespace amrex;

//
// Set the BCs for all the variables EXCEPT pressure or velocity.
//
void
incflo::incflo_set_tracer_bcs (Real time,
                               Vector< std::unique_ptr<MultiFab> > & tracer_in)
{
  BL_PROFILE("incflo::incflo_set_tracer_bcs()");

  if (advect_tracer)
  {
     for (int lev = 0; lev <= finest_level; lev++)
     {
        Box domain(geom[lev].Domain());

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(*tracer_in[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            set_tracer_bcs(time, lev, (*tracer_in[lev])[mfi], 0, tracer_in[lev]->nComp(), domain);
        }

        tracer_in[lev] -> FillBoundary (geom[lev].periodicity());
#ifdef AMREX_USE_EB
        EB_set_covered(*tracer_in[lev], 0, tracer_in[lev]->nComp(), tracer_in[lev]->nGrow(), covered_val);
#endif
     }
  }
}

void 
incflo::set_tracer_bcs(Real time,
                       const int lev,
                       FArrayBox& scal_fab,
                       int bc_comp,
                       int ncomp,
                       const Box& domain)

{
}
