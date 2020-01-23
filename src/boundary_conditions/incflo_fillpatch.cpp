#include <incflo.H>
#include <prob_bc.H>
#include <AMReX_FillPatchUtil.H>
#ifdef AMREX_USE_EB
#include <AMReX_EBInterpolater.H>
#endif

using namespace amrex;

void incflo::fillpatch_velocity (int lev, Real time, MultiFab& vel, int ng)
{
    if (lev == 0) {
        PhysBCFunct<GpuBndryFuncFab<IncfloVelFill> > physbc
            (geom[lev], get_velocity_bcrec(),
             IncfloVelFill{m_probtype, m_bc_velocity});
        FillPatchSingleLevel(vel, IntVect(ng), time,
                             {&(m_leveldata[lev]->velocity_o),
                              &(m_leveldata[lev]->velocity)},
                             {m_t_old[lev], m_t_new[lev]}, 0, 0, 3, geom[lev],
                             physbc, 0);
    } else {
        const auto& bcrec = get_velocity_bcrec();
        PhysBCFunct<GpuBndryFuncFab<IncfloVelFill> > cphysbc
            (geom[lev-1], bcrec, IncfloVelFill{m_probtype, m_bc_velocity});
        PhysBCFunct<GpuBndryFuncFab<IncfloVelFill> > fphysbc
            (geom[lev], bcrec, IncfloVelFill{m_probtype, m_bc_velocity});
#ifdef AMREX_USE_EB
	Interpolater* mapper = (EBFactory(0).isAllRegular()) ?
            (Interpolater*)(&cell_cons_interp) : (Interpolater*)(&eb_cell_cons_interp);
#else
	Interpolater* mapper = &cell_cons_interp;
#endif
        FillPatchTwoLevels(vel, IntVect(ng), time,
                           {&(m_leveldata[lev-1]->velocity_o),
                            &(m_leveldata[lev-1]->velocity)},
                           {m_t_old[lev-1], m_t_new[lev-1]},
                           {&(m_leveldata[lev]->velocity_o),
                            &(m_leveldata[lev]->velocity)},
                           {m_t_old[lev], m_t_new[lev]},
                           0, 0, 3, geom[lev-1], geom[lev],
                           cphysbc, 0, fphysbc, 0,
                           refRatio(lev-1), mapper, bcrec, 0);
    }
}

void incflo::fillpatch_density (int lev, Real time, MultiFab& density, int ng)
{
    if (lev == 0) {
        PhysBCFunct<GpuBndryFuncFab<IncfloDenFill> > physbc(geom[lev], get_density_bcrec(),
                                                            IncfloDenFill{m_probtype, m_bc_density});
        FillPatchSingleLevel(density, IntVect(ng), time,
                             {&(m_leveldata[lev]->density_o),
                              &(m_leveldata[lev]->density)},
                             {m_t_old[lev], m_t_new[lev]}, 0, 0, 1, geom[lev],
                             physbc, 0);
    } else {
        const auto& bcrec = get_density_bcrec();
        PhysBCFunct<GpuBndryFuncFab<IncfloDenFill> > cphysbc
            (geom[lev-1], bcrec, IncfloDenFill{m_probtype, m_bc_density});
        PhysBCFunct<GpuBndryFuncFab<IncfloDenFill> > fphysbc
            (geom[lev], bcrec, IncfloDenFill{m_probtype, m_bc_density});
#ifdef AMREX_USE_EB
	Interpolater* mapper = (EBFactory(0).isAllRegular()) ?
            (Interpolater*)(&cell_cons_interp) : (Interpolater*)(&eb_cell_cons_interp);
#else
	Interpolater* mapper = &cell_cons_interp;
#endif
        FillPatchTwoLevels(density, IntVect(ng), time,
                           {&(m_leveldata[lev-1]->density_o),
                            &(m_leveldata[lev-1]->density)},
                           {m_t_old[lev-1], m_t_new[lev-1]},
                           {&(m_leveldata[lev]->density_o),
                            &(m_leveldata[lev]->density)},
                           {m_t_old[lev], m_t_new[lev]},
                           0, 0, 1, geom[lev-1], geom[lev],
                           cphysbc, 0, fphysbc, 0,
                           refRatio(lev-1), mapper, bcrec, 0);
    }
}

void incflo::fillpatch_tracer (int lev, Real time, MultiFab& tracer, int ng)
{
    if (m_ntrac <= 0) return;
    if (lev == 0) {
        PhysBCFunct<GpuBndryFuncFab<IncfloTracFill> > physbc
            (geom[lev], get_tracer_bcrec(), IncfloTracFill{m_probtype, m_ntrac, m_bc_tracer_d});
        FillPatchSingleLevel(tracer, IntVect(ng), time,
                             {&(m_leveldata[lev]->tracer_o),
                              &(m_leveldata[lev]->tracer)},
                             {m_t_old[lev], m_t_new[lev]}, 0, 0, 1, geom[lev],
                             physbc, 0);
    } else {
        const auto& bcrec = get_tracer_bcrec();
        PhysBCFunct<GpuBndryFuncFab<IncfloTracFill> > cphysbc
            (geom[lev-1], bcrec, IncfloTracFill{m_probtype, m_ntrac, m_bc_tracer_d});
        PhysBCFunct<GpuBndryFuncFab<IncfloTracFill> > fphysbc
            (geom[lev], bcrec, IncfloTracFill{m_probtype, m_ntrac, m_bc_tracer_d});
#ifdef AMREX_USE_EB
	Interpolater* mapper = (EBFactory(0).isAllRegular()) ?
            (Interpolater*)(&cell_cons_interp) : (Interpolater*)(&eb_cell_cons_interp);
#else
	Interpolater* mapper = &cell_cons_interp;
#endif
        FillPatchTwoLevels(tracer, IntVect(ng), time,
                           {&(m_leveldata[lev-1]->tracer_o),
                            &(m_leveldata[lev-1]->tracer)},
                           {m_t_old[lev-1], m_t_new[lev-1]},
                           {&(m_leveldata[lev]->tracer_o),
                            &(m_leveldata[lev]->tracer)},
                           {m_t_old[lev], m_t_new[lev]},
                           0, 0, m_ntrac, geom[lev-1], geom[lev],
                           cphysbc, 0, fphysbc, 0,
                           refRatio(lev-1), mapper, bcrec, 0);
    }
}
