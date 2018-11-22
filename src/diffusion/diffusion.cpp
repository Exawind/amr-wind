#include <AMReX_BC_TYPES.H>
#include <AMReX_Box.H>
#include <AMReX_MLMG.H>
#include <AMReX_MLEBABecLap.H>
#include <AMReX_VisMF.H>

#include <diffusion_F.H>
#include <incflo.H>

//
// Divergence of stress tensor
//
void incflo::ComputeDivTau(int lev,
                           MultiFab& divtau,
                           Vector<std::unique_ptr<MultiFab>>& vel_in)
{
	BL_PROFILE("incflo::ComputeDivTau");
	Box domain(geom[lev].Domain());

   // Get EB geometric info
   Array< const MultiCutFab*,AMREX_SPACEDIM> areafrac;
   Array< const MultiCutFab*,AMREX_SPACEDIM> facecent;
   const amrex::MultiFab*                    volfrac;
   const amrex::MultiCutFab*                 bndrycent;

   areafrac  =   ebfactory[lev] -> getAreaFrac();
   facecent  =   ebfactory[lev] -> getFaceCent();
   volfrac   = &(ebfactory[lev] -> getVolFrac());
   bndrycent = &(ebfactory[lev] -> getBndryCent());

#ifdef _OPENMP
#pragma omp parallel
#endif
   for (MFIter mfi(*vel_in[lev],true); mfi.isValid(); ++mfi) {

      // Tilebox
      Box bx = mfi.tilebox ();

      // this is to check efficiently if this tile contains any eb stuff
      const EBFArrayBox&  vel_fab = static_cast<EBFArrayBox const&>((*vel_in[lev])[mfi]);
      const EBCellFlagFab&  flags = vel_fab.getEBCellFlagFab();

      if (flags.getType(bx) == FabType::covered)
      {
         divtau[mfi].setVal(1.2345e200, bx, 0, 3);
      }
      else
      {
         if (flags.getType(amrex::grow(bx,nghost)) == FabType::regular)
         {
            compute_divtau(
               BL_TO_FORTRAN_BOX(bx),
               BL_TO_FORTRAN_ANYD(divtau[mfi]),
               BL_TO_FORTRAN_ANYD((*vel_in[lev])[mfi]),
               (*eta[lev])[mfi].dataPtr(),
               BL_TO_FORTRAN_ANYD((*ro[lev])[mfi]),
               domain.loVect (), domain.hiVect (),
               bc_ilo[lev]->dataPtr(), bc_ihi[lev]->dataPtr(),
               bc_jlo[lev]->dataPtr(), bc_jhi[lev]->dataPtr(),
               bc_klo[lev]->dataPtr(), bc_khi[lev]->dataPtr(),
               geom[lev].CellSize(), &nghost, &explicit_diffusion);
         }
         else
         {
            compute_divtau_eb(
               BL_TO_FORTRAN_BOX(bx),
               BL_TO_FORTRAN_ANYD(divtau[mfi]),
               BL_TO_FORTRAN_ANYD((*vel_in[lev])[mfi]),
               (*eta[lev])[mfi].dataPtr(),
               BL_TO_FORTRAN_ANYD((*ro[lev])[mfi]),
               BL_TO_FORTRAN_ANYD(flags),
               BL_TO_FORTRAN_ANYD((*areafrac[0])[mfi]),
               BL_TO_FORTRAN_ANYD((*areafrac[1])[mfi]),
               BL_TO_FORTRAN_ANYD((*areafrac[2])[mfi]),
               BL_TO_FORTRAN_ANYD((*facecent[0])[mfi]),
               BL_TO_FORTRAN_ANYD((*facecent[1])[mfi]),
               BL_TO_FORTRAN_ANYD((*facecent[2])[mfi]),
               BL_TO_FORTRAN_ANYD((*volfrac)[mfi]),
               BL_TO_FORTRAN_ANYD((*bndrycent)[mfi]),
               domain.loVect (), domain.hiVect (),
               bc_ilo[lev]->dataPtr(), bc_ihi[lev]->dataPtr(),
               bc_jlo[lev]->dataPtr(), bc_jhi[lev]->dataPtr(),
               bc_klo[lev]->dataPtr(), bc_khi[lev]->dataPtr(),
               geom[lev].CellSize(), &nghost, &explicit_diffusion);

         }
      }
   }
}

//
// Implicit diffusion: solve 
//
//      ( 1 - del dot ( eta grad ) ) u_new = RHS
//
// Here RHS = "vel" which is the current approximation to 
// the new-time velocity (without the diffusion terms solved for implicitly)
//
void incflo::DiffuseVelocity(amrex::Real time)
{
	BL_PROFILE("incflo::DiffuseVelocity");

	// Swap ghost cells and apply BCs to velocity
	FillVelocityBC(time, 0);

	// Update the solver with the relevant current states of the simulation
	diffusion_equation->setCurrentState(ro, eta, dt);

	// Loop over the velocity components
	for(int dir = 0; dir < 3; dir++)
	{
		diffusion_equation->solve(vel, ro, dir);
	}

	// Swap ghost cells and apply BCs to velocity
	FillVelocityBC(time, 0);
}

