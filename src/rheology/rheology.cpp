#include <AMReX_Box.H>

#include <incflo.H>
#include <rheology_F.H>

void incflo::ComputeViscosity()
{
	BL_PROFILE("incflo::ComputeViscosity");

    for(int lev = 0; lev <= finest_level; lev++)
    {
        Box domain(geom[lev].Domain());

#ifdef _OPENMP
#pragma omp parallel
#endif
        for(MFIter mfi(*vel[lev], true); mfi.isValid(); ++mfi)
        {
            // Tilebox
            Box bx = mfi.tilebox();

            const auto& strainrate_arr = strainrate[lev]->array(mfi);
            const auto& viscosity_arr = eta[lev]->array(mfi);

            for(int i = bx.smallEnd(0); i<= bx.bigEnd(0); i++)
            for(int j = bx.smallEnd(1); j<= bx.bigEnd(1); j++)
#if (AMREX_SPACEDIM == 2)
            {
                viscosity_arr(i,j) = viscosity(strainrate_arr(i,j));
            }
#elif (AMREX_SPACEDIM == 3)
            for(int k = bx.smallEnd(2); k<= bx.bigEnd(2); k++)
            {
                viscosity_arr(i,j,k) = viscosity(strainrate_arr(i,j,k));
            }
#endif
        }
    }
}
