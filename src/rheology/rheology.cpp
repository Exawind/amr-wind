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

            if(fluid_model == "newtonian")
            {
                // Viscosity is constant, has already been set in init_fluid()
            }
            else if(fluid_model == "powerlaw")
            {
                powerlaw_viscosity(BL_TO_FORTRAN_BOX(bx),
                                   BL_TO_FORTRAN_ANYD((*eta[lev])[mfi]),
                                   BL_TO_FORTRAN_ANYD((*strainrate[lev])[mfi]));
            }
            else if(fluid_model == "bingham")
            {
                bingham_viscosity(BL_TO_FORTRAN_BOX(bx),
                                  BL_TO_FORTRAN_ANYD((*eta[lev])[mfi]),
                                  BL_TO_FORTRAN_ANYD((*strainrate[lev])[mfi]));
            }
            else if(fluid_model == "hb")
            {
                hb_viscosity(BL_TO_FORTRAN_BOX(bx),
                             BL_TO_FORTRAN_ANYD((*eta[lev])[mfi]),
                             BL_TO_FORTRAN_ANYD((*strainrate[lev])[mfi]));
            }
            else if(fluid_model == "smd")
            {
                smd_viscosity(BL_TO_FORTRAN_BOX(bx),
                              BL_TO_FORTRAN_ANYD((*eta[lev])[mfi]),
                              BL_TO_FORTRAN_ANYD((*strainrate[lev])[mfi]));
            }
            else
            {
                // This should have been caught earlier, but doesn't hurt to double check
                amrex::Abort("Unknown fluid_model! Choose either newtonian, powerlaw, bingham, hb, smd");
            }
        }
    }
}
