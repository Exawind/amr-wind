#include <incflo.H>
#include <prob_bc.H>
#include <AMReX_FillPatchUtil.H>

using namespace amrex;

void incflo::fillpatch_force (Real time, Vector<MultiFab*> const& force, int ng)
{
    const int ncomp = force[0]->nComp();
    const auto& bcrec = get_force_bcrec();
    int lev = 0;
    {
        PhysBCFunct<GpuBndryFuncFab<IncfloForFill> > physbc
            (geom[lev], bcrec, IncfloForFill{m_probtype});
        FillPatchSingleLevel(*force[lev], IntVect(ng), time,
                             {force[lev]}, {time},
                             0, 0, ncomp, geom[lev],
                             physbc, 0);
    }
    for (lev = 1; lev <= finest_level; ++lev)
    {
        PhysBCFunct<GpuBndryFuncFab<IncfloForFill> > cphysbc
            (geom[lev-1], bcrec, IncfloForFill{m_probtype});
        PhysBCFunct<GpuBndryFuncFab<IncfloForFill> > fphysbc
            (geom[lev  ], bcrec, IncfloForFill{m_probtype});
        Interpolater* mapper = &pc_interp;
        FillPatchTwoLevels(*force[lev], IntVect(ng), time,
                           {force[lev-1]}, {time},
                           {force[lev  ]}, {time},
                           0, 0, ncomp, geom[lev-1], geom[lev],
                           cphysbc, 0, fphysbc, 0,
                           refRatio(lev-1), mapper, bcrec, 0);
    }
}
