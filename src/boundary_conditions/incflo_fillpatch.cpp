#include <incflo.H>
#include <prob_bc.H>
// #include <incflo_F.H>
// #include <incflo_util_F.H>
#include <AMReX_FillPatchUtil.H>

using namespace amrex;

void incflo::fillpatch_velocity (int lev, Real time, MultiFab& vel, int ng)
{
    if (lev == 0) {
        PhysBCFunct<GpuBndryFuncFab<IncfloVelFill> > physbc(geom[lev], get_velocity_bcrec(),
                                                            IncfloVelFill{probtype, m_bc_velocity});
        FillPatchSingleLevel(vel, IntVect(ng), time,
                             {&(m_leveldata[lev]->velocity_o),
                              &(m_leveldata[lev]->velocity)},
                             {t_old[lev], t_new[lev]}, 0, 0, 3, geom[lev],
                             physbc, 0);
    } else {
        amrex::Abort("fillpatch_velocity: multi-level todo");
    }
}

void incflo::fillpatch_density (int lev, Real time, MultiFab& density, int ng)
{
    if (lev == 0) {
        PhysBCFunct<GpuBndryFuncFab<IncfloDenFill> > physbc(geom[lev], get_density_bcrec(),
                                                            IncfloDenFill{probtype, m_bc_density});
        FillPatchSingleLevel(density, IntVect(ng), time,
                             {&(m_leveldata[lev]->density_o),
                              &(m_leveldata[lev]->density)},
                             {t_old[lev], t_new[lev]}, 0, 0, 1, geom[lev],
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
            (geom[lev], get_tracer_bcrec(), IncfloTracFill{probtype, ntrac, m_bc_tracer_d});
        FillPatchSingleLevel(tracer, IntVect(ng), time,
                             {&(m_leveldata[lev]->tracer_o),
                              &(m_leveldata[lev]->tracer)},
                             {t_old[lev], t_new[lev]}, 0, 0, 1, geom[lev],
                             physbc, 0);
    } else {
        amrex::Abort("fillpatch_tracer: multi-level todo");
    }
}

namespace
{
    incflo* my_incflo = nullptr;
}

// This interface must match the definition of the interface for
//    CpuBndryFuncFab in amrex/Src/Base/AMReX_PhysBCFunct.H
static
void VelFillBox (Box const& bx, Array4<Real> const& dest,
                 const int dcomp, const int ncomp,
                 GeometryData const& geom, const Real time_in,
                 const BCRec* bcr, const int bcomp,
                 const int orig_comp)
{
    if (dcomp != 0)
         amrex::Abort("Must have dcomp = 0 in VelFillBox");
    if (ncomp != 3)
         amrex::Abort("Must have ncomp = 3 in VelFillBox");

    const Box& domain = geom.Domain();
    int lev = my_incflo->GetLevel(domain);

    // We only do this to make it not const
    Real time = time_in;

    FArrayBox dest_fab(dest);
    Elixir eli_dest_fab = dest_fab.elixir();

    my_incflo->set_velocity_bcs (time, lev, dest_fab, domain);
}

// This interface must match the definition of the interface for
//    CpuBndryFuncFab in amrex/Src/Base/AMReX_PhysBCFunct.H
static
void DensityFillBox (Box const& bx, Array4<Real> const& dest,
                     const int dcomp, const int ncomp,
                     GeometryData const& geom, const Real time_in,
                     const BCRec* bcr, const int bcomp,
                     const int orig_comp)
{
    if (dcomp != 0)
         amrex::Abort("Must have dcomp = 0 in DensityFillBox");
    if (ncomp != 1)
         amrex::Abort("Must have ncomp = 1 in DensityFillBox");

    const Box& domain = geom.Domain();
    int lev = my_incflo->GetLevel(domain);

    // We only do this to make it not const
    Real time = time_in;

    FArrayBox dest_fab(dest);
    Elixir eli_dest_fab = dest_fab.elixir();

    my_incflo->set_density_bcs (time, lev, dest_fab, domain);
}

// This interface must match the definition of the interface for
//    CpuBndryFuncFab in amrex/Src/Base/AMReX_PhysBCFunct.H
static
void ScalarFillBox (Box const& bx, Array4<Real> const& dest,
                    const int dcomp, const int ncomp,
                    GeometryData const& geom, const Real time_in,
                    const BCRec* bcr, const int bcomp,
                    const int orig_comp)
{
    if (dcomp != 0)
         amrex::Abort("Must have dcomp = 0 in ScalaryFillBox");

    const Box& domain = geom.Domain();
    int lev = my_incflo->GetLevel(domain);

    // We only do this to make it not const
    Real time = time_in;

    FArrayBox dest_fab(dest);
    Elixir eli_dest_fab = dest_fab.elixir();

    my_incflo->set_tracer_bcs (time, lev, dest_fab, dcomp, ncomp, domain);
}

// Compute a new multifab by copying array from valid region and filling ghost cells
// works for single level and 2-level cases (fill fine grid ghost by interpolating from coarse)
void
incflo::FillPatchVel (int lev, Real time, MultiFab& mf)
{
    if (my_incflo == nullptr) my_incflo = this;

    // Hack so that ghost cells are not undefined
    mf.setVal(covered_val);

    // There aren't used for anything but need to be defined for the function call
    Vector<BCRec> bcs(3);

    if (lev == 0)
    {
        Vector<MultiFab*> smf;
        Vector<Real> stime;
        GetDataVel(0, time, smf, stime);

        CpuBndryFuncFab bfunc(VelFillBox);
        PhysBCFunct<CpuBndryFuncFab> physbc(geom[lev], bcs, bfunc);
        amrex::FillPatchSingleLevel(mf, time, smf, stime, 0, 0, 3, geom[lev], physbc, 0);
    }
    else
    {
        Vector<MultiFab*> cmf, fmf;
        Vector<Real> ctime, ftime;
        GetDataVel(lev-1, time, cmf, ctime);
        GetDataVel(lev  , time, fmf, ftime);

        CpuBndryFuncFab bfunc(VelFillBox);
        PhysBCFunct<CpuBndryFuncFab> cphysbc(geom[lev-1],bcs,bfunc);
        PhysBCFunct<CpuBndryFuncFab> fphysbc(geom[lev  ],bcs,bfunc);

        Interpolater* mapper = &cell_cons_interp;

        amrex::FillPatchTwoLevels(mf, time, cmf, ctime, fmf, ftime,
                                  0, 0, 3, geom[lev-1], geom[lev],
                                  cphysbc, 0, fphysbc, 0,
                                  refRatio(lev-1), mapper, bcs, 0);
    }
}

// Compute a new multifab by copying array from valid region and filling ghost cells
// works for single level and 2-level cases (fill fine grid ghost by interpolating from coarse)
void
incflo::FillPatchDensity (int lev, Real time, MultiFab& mf)
{
    if (my_incflo == nullptr) my_incflo = this;

    // Hack so that ghost cells are not undefined
    mf.setVal(covered_val);

    // There aren't used for anything but need to be defined for the function call
    Vector<BCRec> bcs(3);

    if (lev == 0)
    {
        Vector<MultiFab*> smf;
        Vector<Real> stime;

        GetDataDensity(0, time, smf, stime);

        CpuBndryFuncFab bfunc(DensityFillBox);
        PhysBCFunct<CpuBndryFuncFab> physbc(geom[lev], bcs, bfunc);
        amrex::FillPatchSingleLevel(mf, time, smf, stime, 0, 0, 1, geom[lev], physbc, 0);
    }
    else
    {
        Vector<MultiFab*> cmf, fmf;
        Vector<Real> ctime, ftime;
        GetDataDensity(lev-1, time, cmf, ctime);
        GetDataDensity(lev  , time, fmf, ftime);

        CpuBndryFuncFab bfunc(DensityFillBox);
        PhysBCFunct<CpuBndryFuncFab> cphysbc(geom[lev-1],bcs,bfunc);
        PhysBCFunct<CpuBndryFuncFab> fphysbc(geom[lev  ],bcs,bfunc);

        Interpolater* mapper = &cell_cons_interp;

        amrex::FillPatchTwoLevels(mf, time, cmf, ctime, fmf, ftime,
                                  0, 0, 1, geom[lev-1], geom[lev],
                                  cphysbc, 0, fphysbc, 0,
                                  refRatio(lev-1), mapper, bcs, 0);
    }
}

// Compute a new multifab by copying array from valid region and filling ghost cells
// works for single level and 2-level cases (fill fine grid ghost by interpolating from coarse)
// NOTE: icomp here refers to whether we are filling 0: density, 1-ntrac-1: tracer(s)
void
incflo::FillPatchScalar (int lev, Real time, MultiFab& mf) 
{
    if (my_incflo == nullptr) my_incflo = this;

    // Hack so that ghost cells are not undefined
    mf.setVal(covered_val);

    // There aren't used for anything but need to be defined for the function call
    Vector<BCRec> bcs(3);

    if (lev == 0)
    {
        Vector<MultiFab*> smf;
        Vector<Real> stime;

        GetDataScalar(0, time, smf, stime);

        CpuBndryFuncFab bfunc(ScalarFillBox);
        PhysBCFunct<CpuBndryFuncFab> physbc(geom[lev], bcs, bfunc);
        amrex::FillPatchSingleLevel(mf, time, smf, stime, 0, 0, ntrac, 
                                    geom[lev], physbc, 0);
    }
    else
    {
        Vector<MultiFab*> cmf, fmf;
        Vector<Real> ctime, ftime;
        GetDataScalar(lev-1, time, cmf, ctime);
        GetDataScalar(lev  , time, fmf, ftime);

        CpuBndryFuncFab bfunc(ScalarFillBox);
        PhysBCFunct<CpuBndryFuncFab> cphysbc(geom[lev-1],bcs,bfunc);
        PhysBCFunct<CpuBndryFuncFab> fphysbc(geom[lev  ],bcs,bfunc);

        Interpolater* mapper = &cell_cons_interp;

        amrex::FillPatchTwoLevels(mf, time, cmf, ctime, fmf, ftime,
                                  0, 0, ntrac, geom[lev-1], geom[lev],
                                  cphysbc, 0, fphysbc, 0,
                                  refRatio(lev-1), mapper, bcs, 0);
    }
}

// Utility to copy in data from phi_old and/or phi_new into another multifab
void
incflo::GetDataVel (int lev, Real time, Vector<MultiFab*>& data, Vector<Real>& datatime)
{
    data.clear();
    datatime.clear();

    const Real teps = (t_new[lev] - t_old[lev]) * 1.e-3;

    if (time > t_new[lev] - teps && time < t_new[lev] + teps)
    {
        data.push_back(vel[lev].get());
        datatime.push_back(t_new[lev]);
    }
    else if (time > t_old[lev] - teps && time < t_old[lev] + teps)
    {
        data.push_back(vel_o[lev].get());
        datatime.push_back(t_old[lev]);
    }
    else
    {
        data.push_back(vel_o[lev].get());
        data.push_back(vel[lev].get());
        datatime.push_back(t_old[lev]);
        datatime.push_back(t_new[lev]);
    }
}

// Utility to copy in data from phi_old and/or phi_new into another multifab
void
incflo::GetDataDensity (int lev, Real time, Vector<MultiFab*>& data, Vector<Real>& datatime)
{
    data.clear();
    datatime.clear();

    const Real teps = (t_new[lev] - t_old[lev]) * 1.e-3;

    if (time > t_new[lev] - teps && time < t_new[lev] + teps)
    {
        data.push_back(density[lev].get());
        datatime.push_back(t_new[lev]);
    }
    else if (time > t_old[lev] - teps && time < t_old[lev] + teps)
    {
        data.push_back(density_o[lev].get());
        datatime.push_back(t_old[lev]);
    }
    else
    {
        data.push_back(density_o[lev].get());
        data.push_back(  density[lev].get());
        datatime.push_back(t_old[lev]);
        datatime.push_back(t_new[lev]);
    }
}

// Utility to copy in data from phi_old and/or phi_new into another multifab
void
incflo::GetDataScalar (int lev, Real time, Vector<MultiFab*>& data, Vector<Real>& datatime)
{
    data.clear();
    datatime.clear();

    const Real teps = (t_new[lev] - t_old[lev]) * 1.e-3;

    if (time > t_new[lev] - teps && time < t_new[lev] + teps)
    {
        data.push_back(tracer[lev].get());
        datatime.push_back(t_new[lev]);
    }
    else if (time > t_old[lev] - teps && time < t_old[lev] + teps)
    {
        data.push_back(tracer_o[lev].get());
        datatime.push_back(t_old[lev]);
    }
    else
    {
        data.push_back(tracer_o[lev].get());
        data.push_back(  tracer[lev].get());
        datatime.push_back(t_old[lev]);
        datatime.push_back(t_new[lev]);
    }
}

