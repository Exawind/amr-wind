#include <AMReX_Box.H>

#include <incflo.H>
#include <rheology_F.H>
#include <boundary_conditions_F.H>

void incflo::ComputeViscosity()
{
    BL_PROFILE("incflo::ComputeViscosity");

    if (fluid_model == "newtonian")
    {
       for(int lev = 0; lev <= finest_level; lev++)
          eta[lev]->setVal(mu,0,1,eta[lev]->nGrow());

    } else {

      // Only compute strain rate if we're going to use it
      ComputeStrainrate();

      for(int lev = 0; lev <= finest_level; lev++)
      {
        Box domain(geom[lev].Domain());

#ifdef _OPENMP
#pragma omp parallel
#endif
        for(MFIter mfi(*eta[lev], true); mfi.isValid(); ++mfi)
        {
            // Tilebox
            Box bx = mfi.tilebox();

            const auto& strainrate_arr = strainrate[lev]->array(mfi);
            const auto& viscosity_arr = eta[lev]->array(mfi);

            AMREX_FOR_3D(bx, i, j, k, 
            {
                viscosity_arr(i,j,k) = viscosity(strainrate_arr(i,j,k));
            });
        }

        eta[lev]->FillBoundary(geom[lev].periodicity());

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for(MFIter mfi(*density[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            fill_bc0(BL_TO_FORTRAN_ANYD((*eta[lev])[mfi]),
                     bc_ilo[lev]->dataPtr(), bc_ihi[lev]->dataPtr(),
                     bc_jlo[lev]->dataPtr(), bc_jhi[lev]->dataPtr(),
                     bc_klo[lev]->dataPtr(), bc_khi[lev]->dataPtr(),
                     domain.loVect(), domain.hiVect(),
                     &nghost);
        }

        eta[lev]->FillBoundary(geom[lev].periodicity());
      }
    }
}
