#include <AMReX_Box.H>
#include <AMReX_EBMultiFabUtil.H>

#include <incflo.H>
#include <derive_F.H>

void incflo::incflo_compute_strainrate(int lev)
{
	BL_PROFILE("incflo::incflo_compute_strainrate");
	Box domain(geom[lev].Domain());

#ifdef _OPENMP
#pragma omp parallel
#endif
	for(MFIter mfi(*vel[lev], true); mfi.isValid(); ++mfi)
	{
		// Tilebox
		Box bx = mfi.tilebox();

		// This is to check efficiently if this tile contains any eb stuff
		const EBFArrayBox& vel_fab = dynamic_cast<EBFArrayBox const&>((*vel[lev])[mfi]);
		const EBCellFlagFab& flags = vel_fab.getEBCellFlagFab();

		if(flags.getType(amrex::grow(bx, 0)) == FabType::regular)
		{
			compute_strainrate(BL_TO_FORTRAN_BOX(bx),
						       BL_TO_FORTRAN_ANYD((*strainrate[lev])[mfi]),
						       BL_TO_FORTRAN_ANYD((*vel[lev])[mfi]),
						       geom[lev].CellSize());
		}
		else
		{
            // TODO: How to compute this in cut cells? 
            // compute_strainrate_eb();
			strainrate[lev]->setVal(0.0, bx, 0, 1);
		}
	}
}

void incflo::incflo_compute_vort(int lev)
{
	BL_PROFILE("incflo::incflo_compute_vort");
	Box domain(geom[lev].Domain());

#ifdef _OPENMP
#pragma omp parallel
#endif
	for(MFIter mfi(*vel[lev], true); mfi.isValid(); ++mfi)
	{
		// Tilebox
		Box bx = mfi.tilebox();

		// This is to check efficiently if this tile contains any eb stuff
		const EBFArrayBox& vel_fab = dynamic_cast<EBFArrayBox const&>((*vel[lev])[mfi]);
		const EBCellFlagFab& flags = vel_fab.getEBCellFlagFab();

		if(flags.getType(amrex::grow(bx, 0)) == FabType::regular)
		{
			compute_vort(BL_TO_FORTRAN_BOX(bx),
						 BL_TO_FORTRAN_ANYD((*vort[lev])[mfi]),
						 BL_TO_FORTRAN_ANYD((*vel[lev])[mfi]),
						 geom[lev].CellSize());
		}
		else
		{
            // TODO: How to compute this in cut cells? 
            // compute_vort_eb();
			vort[lev]->setVal(0.0, bx, 0, 1);
		}
	}
}
