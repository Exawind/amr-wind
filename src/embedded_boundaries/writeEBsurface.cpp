#include <incflo.H>
#include <AMReX_WriteEBSurface.H>

using namespace amrex;

void incflo::WriteMyEBSurface()
{
  if (Geom(0).isAllPeriodic()) return;

  // Only write at the finest level!
  int lev = 0;

  BoxArray & ba            = grids[lev];
  DistributionMapping & dm = dmap[lev];

  const EBFArrayBoxFactory * ebf;

  ebf = ebfactory[lev].get();

  WriteEBSurface(ba,dm,Geom(lev),ebf);
}
