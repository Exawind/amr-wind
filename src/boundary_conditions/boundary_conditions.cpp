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
incflo::FillPatch (int lev, MultiFab& mf, MultiFab& cmf, MultiFab& fmf, int icomp, int ncomp)
{
#if 0
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
#endif
}

//
// Set the BCs for all the variables EXCEPT pressure or velocity.
//
void incflo::incflo_set_scalar_bcs ()
{
    BL_PROFILE("incflo::incflo_set_scalar_bcs()");

    for(int lev = 0; lev < nlev; lev++)
    {
        Box domain(geom[lev].Domain());

#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(*ro[lev], true); mfi.isValid(); ++mfi)
        {
            set_scalar_bcs((*ro[lev])[mfi].dataPtr(),
                           BL_TO_FORTRAN_ANYD((*eta[lev])[mfi]),
                           bc_ilo[lev]->dataPtr(), bc_ihi[lev]->dataPtr(),
                           bc_jlo[lev]->dataPtr(), bc_jhi[lev]->dataPtr(),
                           bc_klo[lev]->dataPtr(), bc_khi[lev]->dataPtr(),
                           domain.loVect(), domain.hiVect(),
                           &nghost);
        }
        ro[lev]->FillBoundary(geom[lev].periodicity());
        eta[lev]->FillBoundary(geom[lev].periodicity());
    }
}

//
// Set the BCs for velocity only
//
void incflo::incflo_set_velocity_bcs(Real time, int extrap_dir_bcs)
{
    BL_PROFILE("incflo::incflo_set_velocity_bcs()");

    for(int lev = 0; lev < nlev; lev++)
    {
        vel[lev]->FillBoundary(geom[lev].periodicity());
        Box domain(geom[lev].Domain());
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
    }
}

//
// Fills ghost cell values of pressure appropriately for the BC type
//
void incflo::incflo_extrap_pressure (int lev, std::unique_ptr<amrex::MultiFab>& p)
{
    BL_PROFILE("incflo::incflo_extrap_pressure()");
 
    Box domain(geom[lev].Domain());
 
#ifdef _OPENMP
#pragma omp parallel
#endif
    for(MFIter mfi(*p, true); mfi.isValid(); ++mfi) 
    {
        extrap_pressure_to_ghost_cells(BL_TO_FORTRAN_ANYD((*p)[mfi]),
                                       bc_ilo[lev]->dataPtr(), bc_ihi[lev]->dataPtr(),
                                       bc_jlo[lev]->dataPtr(), bc_jhi[lev]->dataPtr(),
                                       bc_klo[lev]->dataPtr(), bc_khi[lev]->dataPtr(),
                                       domain.loVect(), domain.hiVect(),
                                       &nghost);
    }
}

void incflo::fill_mf_bc(int lev, MultiFab& mf)
{
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

