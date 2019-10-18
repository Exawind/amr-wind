#include <incflo.H>
// #include <incflo_F.H>
// #include <incflo_util_F.H>
#include <AMReX_FillPatchUtil.H>

namespace
{
  incflo* incflo_for_fillpatching;
}

// This interface must match the definition of the interface for
//    CpuBndryFuncFab in amrex/Src/Base/AMReX_PhysBCFunct.H
void set_ptr_to_incflo(incflo& incflo_for_fillpatching_in)
{
   incflo_for_fillpatching = &incflo_for_fillpatching_in;
}

// This interface must match the definition of the interface for
//    CpuBndryFuncFab in amrex/Src/Base/AMReX_PhysBCFunct.H
inline
void VelFillBox (Box const& bx, Array4<amrex::Real> const& dest,
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

    // This is a bit hack-y but does get us the right level
    int lev = 0;
    for (int ilev = 0; ilev < 10; ilev++)
    {
       const Geometry& lev_geom = incflo_for_fillpatching->get_geom_ref(ilev);
       if (domain.length()[0] == (lev_geom.Domain()).length()[0])
       {
         lev = ilev;
         break;
       }
    }

    // We are hard-wiring this fillpatch routine to define the Dirichlet values
    //    at the faces (not the ghost cell center)
    int extrap_dir_bcs = 0;

    // We only do this to make it not const
    Real time = time_in;

    FArrayBox dest_fab(dest);

    incflo_for_fillpatching->set_velocity_bcs (time, lev, dest_fab, domain, &extrap_dir_bcs);
}

// This interface must match the definition of the interface for
//    CpuBndryFuncFab in amrex/Src/Base/AMReX_PhysBCFunct.H
inline
void DensityFillBox (Box const& bx, Array4<amrex::Real> const& dest,
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

    // This is a bit hack-y but does get us the right level
    int lev = 0;
    for (int ilev = 0; ilev < 10; ilev++)
    {
       const Geometry& lev_geom = incflo_for_fillpatching->get_geom_ref(ilev);
       if (domain.length()[0] == (lev_geom.Domain()).length()[0])
       {
         lev = ilev;
         break;
       }
    }

    // We only do this to make it not const
    Real time = time_in;

    FArrayBox dest_fab(dest);

    incflo_for_fillpatching->set_density_bcs (time, lev, dest_fab, domain);
}

// This interface must match the definition of the interface for
//    CpuBndryFuncFab in amrex/Src/Base/AMReX_PhysBCFunct.H
inline
void ScalarFillBox (Box const& bx, Array4<amrex::Real> const& dest,
                    const int dcomp, const int ncomp,
                    GeometryData const& geom, const Real time_in,
                    const BCRec* bcr, const int bcomp,
                    const int orig_comp)
{
    if (dcomp != 0)
         amrex::Abort("Must have dcomp = 0 in ScalaryFillBox");

    const Box& domain = geom.Domain();

    // This is a bit hack-y but does get us the right level
    int lev = 0;
    for (int ilev = 0; ilev < 10; ilev++)
    {
       const Geometry& lev_geom = incflo_for_fillpatching->get_geom_ref(ilev);
       if (domain.length()[0] == (lev_geom.Domain()).length()[0])
       {
         lev = ilev;
         break;
       }
    }

    // We only do this to make it not const
    Real time = time_in;

    FArrayBox dest_fab(dest);

    incflo_for_fillpatching->set_tracer_bcs (time, lev, dest_fab, dcomp, ncomp, domain);
}

// Compute a new multifab by copying array from valid region and filling ghost cells
// works for single level and 2-level cases (fill fine grid ghost by interpolating from coarse)
void
incflo::FillPatchVel (int lev, Real time, MultiFab& mf)
{
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

