#include <AMReX_Box.H>

#include <incflo.H>
#include <rheology_F.H>

void incflo::incflo_compute_viscosity()
{
	BL_PROFILE("incflo::incflo_compute_viscosity");

    for(int lev = 0; lev < nlev; lev++)
    {
        Box domain(geom[lev].Domain());

    #ifdef _OPENMP
    #pragma omp parallel
    #endif
        for(MFIter mfi(*vel[lev], true); mfi.isValid(); ++mfi)
        {
            // Tilebox
            Box bx = mfi.tilebox();

            compute_viscosity(BL_TO_FORTRAN_BOX(bx),
                              BL_TO_FORTRAN_ANYD((*mu[lev])[mfi]),
                              BL_TO_FORTRAN_ANYD((*strainrate[lev])[mfi]),
                              geom[lev].CellSize());
        }
    }
}
