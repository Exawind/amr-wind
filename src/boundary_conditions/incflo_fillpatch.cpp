#include <incflo.H>
#include <prob_bc.H>
#include <AMReX_FillPatchUtil.H>

using namespace amrex;

void incflo::fillpatch_velocity (int lev, Real time, MultiFab& vel, int ng)
{
    if (lev == 0) {
        PhysBCFunct<GpuBndryFuncFab<IncfloVelFill> > physbc(geom[lev], get_velocity_bcrec(),
                                                            IncfloVelFill{m_probtype, m_bc_velocity});
        FillPatchSingleLevel(vel, IntVect(ng), time,
                             {&(m_leveldata[lev]->velocity_o),
                              &(m_leveldata[lev]->velocity)},
                             {m_t_old[lev], m_t_new[lev]}, 0, 0, 3, geom[lev],
                             physbc, 0);
    } else {
        amrex::Abort("fillpatch_velocity: multi-level todo");
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
        amrex::Abort("fillpatch_density: multi-level todo");
    }
}

void incflo::fillpatch_tracer (int lev, Real time, MultiFab& tracer, int ng)
{
    if (ntrac <= 0) return;
    if (lev == 0) {
        PhysBCFunct<GpuBndryFuncFab<IncfloTracFill> > physbc
            (geom[lev], get_tracer_bcrec(), IncfloTracFill{m_probtype, ntrac, m_bc_tracer_d});
        FillPatchSingleLevel(tracer, IntVect(ng), time,
                             {&(m_leveldata[lev]->tracer_o),
                              &(m_leveldata[lev]->tracer)},
                             {m_t_old[lev], m_t_new[lev]}, 0, 0, 1, geom[lev],
                             physbc, 0);
    } else {
        amrex::Abort("fillpatch_tracer: multi-level todo");
    }
}
