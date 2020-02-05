#include <incflo.H>
#include <prob_bc.H>

using namespace amrex;

void incflo::fillphysbc_velocity (int lev, Real time, MultiFab& vel, int ng)
{
    PhysBCFunct<GpuBndryFuncFab<IncfloVelFill> > physbc(geom[lev], get_velocity_bcrec(),
                                                        IncfloVelFill{m_probtype, m_bc_velocity});
    physbc.FillBoundary(vel, 0, AMREX_SPACEDIM, IntVect(ng), time, 0);
}

void incflo::fillphysbc_density (int lev, Real time, MultiFab& density, int ng)
{
    PhysBCFunct<GpuBndryFuncFab<IncfloDenFill> > physbc(geom[lev], get_density_bcrec(),
                                                        IncfloDenFill{m_probtype, m_bc_density});
    physbc.FillBoundary(density, 0, 1, IntVect(ng), time, 0);
}

void incflo::fillphysbc_tracer (int lev, Real time, MultiFab& tracer, int ng)
{
    if (m_ntrac > 0) {
        PhysBCFunct<GpuBndryFuncFab<IncfloTracFill> > physbc
            (geom[lev], get_tracer_bcrec(), IncfloTracFill{m_probtype, m_ntrac, m_bc_tracer_d});
        physbc.FillBoundary(tracer, 0, m_ntrac, IntVect(ng), time, 0);
    }
}
