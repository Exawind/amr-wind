#include <AMReX_Box.H>
#include <AMReX_EBMultiFabUtil.H>

#include <incflo.H>
#include <derive_F.H>

void incflo::incflo_compute_strainrate()
{
    BL_PROFILE("incflo::incflo_compute_strainrate");

    for(int lev = 0; lev < nlev; lev++)
    {
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
        for(MFIter mfi(*vel[lev], true); mfi.isValid(); ++mfi)
        {
            // Tilebox
            Box bx = mfi.tilebox();

            // This is to check efficiently if this tile contains any eb stuff
            const EBFArrayBox& vel_fab = static_cast<EBFArrayBox const&>((*vel[lev])[mfi]);
            const EBCellFlagFab& flags = vel_fab.getEBCellFlagFab();

            if (flags.getType(bx) == FabType::covered)
            {
                (*strainrate[lev])[mfi].setVal(1.2345e200, bx, 0, 3);
            }
            else
            {
                if(flags.getType(amrex::grow(bx, 0)) == FabType::regular)
                {
                    compute_strainrate(BL_TO_FORTRAN_BOX(bx),
                                       BL_TO_FORTRAN_ANYD((*strainrate[lev])[mfi]),
                                       BL_TO_FORTRAN_ANYD((*vel[lev])[mfi]),
                                       geom[lev].CellSize());
                }
                else
                {
                    compute_strainrate_eb(BL_TO_FORTRAN_BOX(bx),
                                          BL_TO_FORTRAN_ANYD((*strainrate[lev])[mfi]),
                                          BL_TO_FORTRAN_ANYD((*vel[lev])[mfi]),
                                          BL_TO_FORTRAN_ANYD(flags),
                                          BL_TO_FORTRAN_ANYD((*areafrac[0])[mfi]),
                                          BL_TO_FORTRAN_ANYD((*areafrac[1])[mfi]),
                                          BL_TO_FORTRAN_ANYD((*areafrac[2])[mfi]),
                                          BL_TO_FORTRAN_ANYD((*facecent[0])[mfi]),
                                          BL_TO_FORTRAN_ANYD((*facecent[1])[mfi]),
                                          BL_TO_FORTRAN_ANYD((*facecent[2])[mfi]),
                                          BL_TO_FORTRAN_ANYD((*volfrac)[mfi]),
                                          BL_TO_FORTRAN_ANYD((*bndrycent)[mfi]),
                                          geom[lev].CellSize());
                }
            }
        }
    }
}

void incflo::incflo_compute_vort()
{
	BL_PROFILE("incflo::incflo_compute_vort");

    for(int lev = 0; lev < nlev; lev++)
    {
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
        for(MFIter mfi(*vel[lev], true); mfi.isValid(); ++mfi)
        {
            // Tilebox
            Box bx = mfi.tilebox();

            // This is to check efficiently if this tile contains any eb stuff
            const EBFArrayBox& vel_fab = static_cast<EBFArrayBox const&>((*vel[lev])[mfi]);
            const EBCellFlagFab& flags = vel_fab.getEBCellFlagFab();

            if (flags.getType(bx) == FabType::covered)
            {
                (*vort[lev])[mfi].setVal(1.2345e200, bx, 0, 3);
            }
            else
            {
                if(flags.getType(amrex::grow(bx, 0)) == FabType::regular)
                {
                    compute_vort(BL_TO_FORTRAN_BOX(bx),
                                 BL_TO_FORTRAN_ANYD((*vort[lev])[mfi]),
                                 BL_TO_FORTRAN_ANYD((*vel[lev])[mfi]),
                                 geom[lev].CellSize());
                }
                else
                {
                    compute_vort_eb(BL_TO_FORTRAN_BOX(bx),
                                    BL_TO_FORTRAN_ANYD((*vort[lev])[mfi]),
                                    BL_TO_FORTRAN_ANYD((*vel[lev])[mfi]),
                                    BL_TO_FORTRAN_ANYD(flags),
                                    BL_TO_FORTRAN_ANYD((*areafrac[0])[mfi]),
                                    BL_TO_FORTRAN_ANYD((*areafrac[1])[mfi]),
                                    BL_TO_FORTRAN_ANYD((*areafrac[2])[mfi]),
                                    BL_TO_FORTRAN_ANYD((*facecent[0])[mfi]),
                                    BL_TO_FORTRAN_ANYD((*facecent[1])[mfi]),
                                    BL_TO_FORTRAN_ANYD((*facecent[2])[mfi]),
                                    BL_TO_FORTRAN_ANYD((*volfrac)[mfi]),
                                    BL_TO_FORTRAN_ANYD((*bndrycent)[mfi]),
                                    geom[lev].CellSize());
                }
            }

        }
    }
}
