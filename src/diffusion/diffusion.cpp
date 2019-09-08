#include <AMReX_BC_TYPES.H>
#include <AMReX_Box.H>
#include <AMReX_EBMultiFabUtil.H>
#include <AMReX_MLMG.H>
#include <AMReX_MLEBABecLap.H>
#include <AMReX_VisMF.H>

#include <diffusion_F.H>
#include <incflo.H>

//
// Explicit part of divergence of stress tensor: 
// div ( eta (grad u)^T )
//
void incflo::ComputeDivTau(int lev,
                           MultiFab& divtau_in,
                           Vector<std::unique_ptr<MultiFab>>& vel_in)
{
    BL_PROFILE("incflo::ComputeDivTau");
    Box domain(geom[lev].Domain());

    EB_set_covered(*vel[lev], covered_val);

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
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*vel_in[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        // Tilebox
        Box bx = mfi.tilebox();

        // this is to check efficiently if this tile contains any eb stuff
        const EBFArrayBox&  vel_fab = static_cast<EBFArrayBox const&>((*vel_in[lev])[mfi]);
        const EBCellFlagFab&  flags = vel_fab.getEBCellFlagFab();

        if (flags.getType(bx) == FabType::covered)
        {
            divtau_in[mfi].setVal(1.2345e200, bx, 0, AMREX_SPACEDIM);
        }
        else
        {
            if (flags.getType(amrex::grow(bx,nghost)) == FabType::regular)
            {
                compute_divtau(BL_TO_FORTRAN_BOX(bx),
                               BL_TO_FORTRAN_ANYD(divtau_in[mfi]),
                               BL_TO_FORTRAN_ANYD((*vel_in[lev])[mfi]),
                               (*eta[lev])[mfi].dataPtr(),
                               BL_TO_FORTRAN_ANYD((*ro[lev])[mfi]),
                               domain.loVect (), domain.hiVect (),
                               bc_ilo[lev]->dataPtr(), bc_ihi[lev]->dataPtr(),
                               bc_jlo[lev]->dataPtr(), bc_jhi[lev]->dataPtr(),
                               bc_klo[lev]->dataPtr(), bc_khi[lev]->dataPtr(),
                               geom[lev].CellSize(), &nghost);
            }
            else
            {
                compute_divtau_eb(BL_TO_FORTRAN_BOX(bx),
                                  BL_TO_FORTRAN_ANYD(divtau_in[mfi]),
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
                                  geom[lev].CellSize(), &nghost, &cyl_speed);
            }
        }
   }
}
