#include <incflo.H>

using namespace amrex;

//
// Set the BCs for all the variables EXCEPT pressure or velocity.
//
void
incflo::incflo_set_density_bcs (Real time,
                                Vector< std::unique_ptr<MultiFab> > & density_in)
{
  BL_PROFILE("incflo::incflo_set_density_bcs()");

  for (int lev = 0; lev <= finest_level; lev++)
  {
     Box domain(geom[lev].Domain());

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
     for (MFIter mfi(*density_in[lev], true); mfi.isValid(); ++mfi) {
         set_density_bcs(time, lev, (*density_in[lev])[mfi], domain);
     }

     density_in[lev] -> FillBoundary (geom[lev].periodicity());
#ifdef AMREX_USE_EB
     EB_set_covered(*density_in[lev], 0, density_in[lev]->nComp(), density_in[lev]->nGrow(), covered_val);
#endif
  }
}

void 
incflo::set_density_bcs(Real time,
                        const int lev,
                        FArrayBox& scal_fab,
                        const Box& domain)

{
}
