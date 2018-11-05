#include <AMReX_Array.H>
#include <AMReX_BC_TYPES.H>
#include <AMReX_BLassert.H>
#include <AMReX_EBMultiFabUtil.H>
#include <AMReX_MLEBABecLap.H>
#include <AMReX_MLMG.H>
#include <AMReX_MLNodeLaplacian.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_VisMF.H>

#include <incflo.H>
#include <boundary_conditions_F.H>
#include <divergence_F.H>
#include <projection_F.H>

//
// Compute div(u)
//
void incflo::incflo_compute_divu(Real time)
{
    if (nodal_pressure == 1)
    {
        int extrap_dir_bcs = 0;
        incflo_set_velocity_bcs (time, extrap_dir_bcs);

        for (int lev = 0; lev < nlev; lev++)
        {
            Box domain(geom[lev].Domain());
            vel[lev]->FillBoundary (geom[lev].periodicity());     
#ifdef _OPENMP
#pragma omp parallel
#endif
           // Extrapolate Dirichlet values to ghost cells -- but do it differently in that 
           //  no-slip walls are treated exactly like slip walls -- this is only relevant
           //  when going into the projection
           for (MFIter mfi((*vel[lev]), true); mfi.isValid(); ++mfi)
           {
                set_vec_bcs ( BL_TO_FORTRAN_ANYD((*vel[lev])[mfi]),
                              bc_ilo[lev]->dataPtr(), bc_ihi[lev]->dataPtr(),
                              bc_jlo[lev]->dataPtr(), bc_jhi[lev]->dataPtr(),
                              bc_klo[lev]->dataPtr(), bc_khi[lev]->dataPtr(),
                              domain.loVect(), domain.hiVect(),
                              &nghost);
           }

           vel[lev]->FillBoundary (geom[lev].periodicity());
       }

        // Define the operator in order to compute the multi-level divergence
        //
        //        (del dot b sigma grad)) phi
        //
        LPInfo                       info;
        MLNodeLaplacian              matrix(geom, grids, dmap, info);

        // Set domain BCs for Poisson's solver
        // The domain BCs refer to level 0 only
        int bc_lo[3], bc_hi[3];
        Box domain(geom[0].Domain());

        set_ppe_bc(bc_lo, bc_hi,
                   domain.loVect(), domain.hiVect(),
                   &nghost,
                   bc_ilo[0]->dataPtr(), bc_ihi[0]->dataPtr(),
                   bc_jlo[0]->dataPtr(), bc_jhi[0]->dataPtr(),
                   bc_klo[0]->dataPtr(), bc_khi[0]->dataPtr());

        matrix.setDomainBC ( {(LinOpBCType)bc_lo[0], (LinOpBCType)bc_lo[1], (LinOpBCType)bc_lo[2]},
                             {(LinOpBCType)bc_hi[0], (LinOpBCType)bc_hi[1], (LinOpBCType)bc_hi[2]} );

        matrix.compDivergence(GetVecOfPtrs(divu), GetVecOfPtrs(vel)); 

    }
    else
    {
        int extrap_dir_bcs = 1;
        incflo_set_velocity_bcs(time, extrap_dir_bcs);

        for(int lev = 0; lev < nlev; lev++)
        {
            Box domain(geom[lev].Domain());
            vel[lev]->FillBoundary(geom[lev].periodicity());

            // Create face centered multifabs for vel
            Array<std::unique_ptr<MultiFab>,AMREX_SPACEDIM> vel_fc;
            incflo_average_cc_to_fc(lev, *vel[lev], vel_fc);

            // This does not need to have correct ghost values in place
            EB_computeDivergence(*divu[lev], GetArrOfConstPtrs(vel_fc), geom[lev]);
        }
    }

	// Restore velocities to carry Dirichlet values on faces
	int extrap_dir_bcs = 0;
	incflo_set_velocity_bcs(time, extrap_dir_bcs);
}

//
// This subroutines averages component by component
// The assumption is that cc is multicomponent
// 
void
incflo::incflo_average_cc_to_fc(int lev, 
                                const MultiFab& cc,
                                Array<std::unique_ptr<MultiFab>,AMREX_SPACEDIM>& fc )
{
    AMREX_ASSERT(cc.nComp()==AMREX_SPACEDIM);
    AMREX_ASSERT(AMREX_SPACEDIM==3);
    
    // First allocate fc
    BoxArray x_ba = cc.boxArray();
    x_ba.surroundingNodes(0);
    fc[0].reset(new MultiFab(x_ba,cc.DistributionMap(),1,nghost, MFInfo(), *ebfactory[lev]));
    fc[0]->setVal(1.0e200);

    BoxArray y_ba = cc.boxArray();
    y_ba.surroundingNodes(1);
    fc[1].reset(new MultiFab(y_ba,cc.DistributionMap(),1,nghost, MFInfo(), *ebfactory[lev]));
    fc[1]->setVal(1.0e200);

    BoxArray z_ba = cc.boxArray();
    z_ba.surroundingNodes(2);
    fc[2].reset(new MultiFab(z_ba,cc.DistributionMap(),1,nghost, MFInfo(), *ebfactory[lev]));
    fc[2]->setVal(1.0e200);

    //
    // Average
    // We do not care about EB because faces in covered regions
    // should never get used so we can set them to whatever values
    // we like
    //
#ifdef _OPENMP
#pragma omp parallel 
#endif
    for (MFIter mfi(*vel[lev],true); mfi.isValid(); ++mfi)
    {
        // Boxes for staggered components
        Box bx = mfi.tilebox();

        average_cc_to_fc( BL_TO_FORTRAN_BOX(bx),
                          BL_TO_FORTRAN_ANYD((*fc[0])[mfi]),
                          BL_TO_FORTRAN_ANYD((*fc[1])[mfi]),
                          BL_TO_FORTRAN_ANYD((*fc[2])[mfi]),
                          BL_TO_FORTRAN_ANYD(cc[mfi]));
    }

    // Set BCs for the fac-centered velocities
	fc[0]->FillBoundary(geom[lev].periodicity());
	fc[1]->FillBoundary(geom[lev].periodicity());
	fc[2]->FillBoundary(geom[lev].periodicity());
    // We do not fill BCs and halo regions in this routine
} 
