#include <AMReX_Array.H>
#include <AMReX_BC_TYPES.H>
#include <AMReX_BLassert.H>
#include <AMReX_Box.H>
#include <AMReX_EBMultiFabUtil.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_MultiFab.H>
#include <AMReX_VisMF.H>

#include <incflo.H>
#include <boundary_conditions_F.H>

// Compute a new multifab by coping in phi from valid region and filling ghost cells
// works for single level and 2-level cases (fill fine grid ghost by interpolating from coarse)
void
incflo::FillPatch (int lev, Real time, MultiFab& mf, MultiFab& cmf, MultiFab& fmf, int icomp, int ncomp)
{
    /*
    if (lev == 0)
    {
        Vector<MultiFab*> smf;
        Vector<Real> stime;
        GetData(0, time, smf, stime);

        PhysBCFunct physbc(geom[lev],bcs,BndryFunctBase(phifill));
        amrex::FillPatchSingleLevel(mf, time, smf, stime, 0, icomp, ncomp,
                                     geom[lev], physbc);
    }
    else
    {
        PhysBCFunct cphysbc(geom[lev-1],bcs,BndryFunctBase(phifill));
        PhysBCFunct fphysbc(geom[lev  ],bcs,BndryFunctBase(phifill));

        Interpolater* mapper = &cell_cons_interp;

        amrex::FillPatchTwoLevels(mf, time, cmf, ctime, fmf, ftime,
                                   0, icomp, ncomp, geom[lev-1], geom[lev],
                                   cphysbc, fphysbc, refRatio(lev-1),
                                   mapper, bcs);
    }
    */
}

//
// Fill the BCs for velocity only
//
void incflo::FillVelocityBC(Real time, int extrap_dir_bcs)
{
    BL_PROFILE("incflo::FillVelocityBC()");

    for(int lev = 0; lev <= finest_level; lev++)
    {
        Box domain(geom[lev].Domain());
        vel[lev]->FillBoundary(geom[lev].periodicity());
#ifdef _OPENMP
#pragma omp parallel
#endif
        for(MFIter mfi(*vel[lev], true); mfi.isValid(); ++mfi)
        {
            set_velocity_bcs(&time, 
                             BL_TO_FORTRAN_ANYD((*vel[lev])[mfi]),
                             bc_ilo[lev]->dataPtr(), bc_ihi[lev]->dataPtr(),
                             bc_jlo[lev]->dataPtr(), bc_jhi[lev]->dataPtr(),
                             bc_klo[lev]->dataPtr(), bc_khi[lev]->dataPtr(),
                             domain.loVect(), domain.hiVect(),
                             &nghost, &extrap_dir_bcs);
        }
        // Do this after as well as before to pick up terms that got updated in the call above
        vel[lev]->FillBoundary(geom[lev].periodicity());
    }
}

void incflo::FillScalarBC(int lev, MultiFab& mf)
{
    BL_PROFILE("incflo:FillScalarBC()");
	Box domain(geom[lev].Domain());

	if(!mf.boxArray().ixType().cellCentered())
		amrex::Error("fill_mf_bc only used for cell-centered arrays!");

	// Impose periodic bc's at domain boundaries and fine-fine copies in the interior
	mf.FillBoundary(geom[lev].periodicity());

    // Fill all cell-centered arrays with first-order extrapolation at domain boundaries
#ifdef _OPENMP
#pragma omp parallel
#endif
	for(MFIter mfi(mf, true); mfi.isValid(); ++mfi)
	{
		const Box& sbx = mf[mfi].box();
		fill_bc0(mf[mfi].dataPtr(),
				 sbx.loVect(),
				 sbx.hiVect(),
				 bc_ilo[lev]->dataPtr(),
				 bc_ihi[lev]->dataPtr(),
				 bc_jlo[lev]->dataPtr(),
				 bc_jhi[lev]->dataPtr(),
				 bc_klo[lev]->dataPtr(),
				 bc_khi[lev]->dataPtr(),
				 domain.loVect(),
				 domain.hiVect(),
				 &nghost);
	}
}

