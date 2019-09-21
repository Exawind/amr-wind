#include <AMReX_Array.H>
#include <AMReX_BC_TYPES.H>
#include <AMReX_BLassert.H>
#include <AMReX_Box.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_VisMF.H>

#include <incflo.H>
#include <boundary_conditions_F.H>
#include <setup_F.H>

namespace
{
  incflo* incflo_for_fillpatching;
}

void set_ptr_to_incflo(incflo& incflo_for_fillpatching_in)
{
   incflo_for_fillpatching = &incflo_for_fillpatching_in;
}

inline void DensityFillBox(Box const& bx, Array4<amrex::Real> const& dest, 
                           const int dcomp, const int numcomp,
                           GeometryData const& geom, const Real time_in, 
                           const BCRec* bcr, 
                           const int bcomp, const int orig_comp)
{
    if (dcomp != 0)
         amrex::Abort("Must have dcomp = 0 in DensityFillBox");
    
    if (numcomp != 1)
         amrex::Abort("Must have numcomp = 3 in DensityFillBox");

    const Box& domain = geom.Domain();

    // This is a bit hack-y but does get us the right level 
    int lev = 0;
    while(lev < 20)
    {
       const Geometry& lev_geom = incflo_for_fillpatching->get_geom_ref(lev);
       if (domain.length()[0] == (lev_geom.Domain()).length()[0]) 
       {
         break;
       }
       lev++;
    }
    if (lev == 20)
        amrex::Abort("Reached lev = 20 in DensityFillBox...");

    // We are hard-wiring this fillpatch routine to define the Dirichlet values
    //    at the faces (not the ghost cell center)
    int extrap_dir_bcs = 1;

    // We only do this to make it not const
    Real time = time_in;

    const int* bc_ilo_ptr = incflo_for_fillpatching->get_bc_ilo_ptr(lev);
    const int* bc_ihi_ptr = incflo_for_fillpatching->get_bc_ihi_ptr(lev);
    const int* bc_jlo_ptr = incflo_for_fillpatching->get_bc_jlo_ptr(lev);
    const int* bc_jhi_ptr = incflo_for_fillpatching->get_bc_jhi_ptr(lev);
    const int* bc_klo_ptr = incflo_for_fillpatching->get_bc_klo_ptr(lev);
    const int* bc_khi_ptr = incflo_for_fillpatching->get_bc_khi_ptr(lev);

    int nghost = incflo_for_fillpatching->get_nghost();
    int probtype = incflo_for_fillpatching->get_probtype();

    FArrayBox dest_fab(dest);

    set_density_bcs(&time, 
                    dest_fab.dataPtr(), dest_fab.loVect(), dest_fab.hiVect(),
                    bc_ilo_ptr, bc_ihi_ptr, 
                    bc_jlo_ptr, bc_jhi_ptr, 
                    bc_klo_ptr, bc_khi_ptr, 
                    domain.loVect(), domain.hiVect(),
                    &nghost, &extrap_dir_bcs, &probtype);
}

inline void TracerFillBox(Box const& bx, Array4<amrex::Real> const& dest, 
                          const int dcomp, const int numcomp,
                          GeometryData const& geom, const Real time_in, 
                          const BCRec* bcr, 
                          const int bcomp, const int orig_comp)
{
    if (dcomp != 0)
         amrex::Abort("Must have dcomp = 0 in TracerFillBox");
    
    if (numcomp != 1)
         amrex::Abort("Must have numcomp = 3 in TracerFillBox");

    const Box& domain = geom.Domain();

    // This is a bit hack-y but does get us the right level 
    int lev = 0;
    while(lev < 20)
    {
       const Geometry& lev_geom = incflo_for_fillpatching->get_geom_ref(lev);
       if (domain.length()[0] == (lev_geom.Domain()).length()[0]) 
       {
         break;
       }
       lev++;
    }
    if (lev == 20)
        amrex::Abort("Reached lev = 20 in TracerFillBox...");

    // We are hard-wiring this fillpatch routine to define the Dirichlet values
    //    at the faces (not the ghost cell center)
    int extrap_dir_bcs = 1;

    // We only do this to make it not const
    Real time = time_in;

    const int* bc_ilo_ptr = incflo_for_fillpatching->get_bc_ilo_ptr(lev);
    const int* bc_ihi_ptr = incflo_for_fillpatching->get_bc_ihi_ptr(lev);
    const int* bc_jlo_ptr = incflo_for_fillpatching->get_bc_jlo_ptr(lev);
    const int* bc_jhi_ptr = incflo_for_fillpatching->get_bc_jhi_ptr(lev);
    const int* bc_klo_ptr = incflo_for_fillpatching->get_bc_klo_ptr(lev);
    const int* bc_khi_ptr = incflo_for_fillpatching->get_bc_khi_ptr(lev);

    int nghost = incflo_for_fillpatching->get_nghost();
    int probtype = incflo_for_fillpatching->get_probtype();

    FArrayBox dest_fab(dest);

    set_tracer_bcs(&time, 
                   dest_fab.dataPtr(), dest_fab.loVect(), dest_fab.hiVect(),
                   bc_ilo_ptr, bc_ihi_ptr, 
                   bc_jlo_ptr, bc_jhi_ptr, 
                   bc_klo_ptr, bc_khi_ptr, 
                   domain.loVect(), domain.hiVect(),
                   &nghost, &extrap_dir_bcs, &probtype);
}

// This interface must match the definition of the interface for
//    CpuBndryFuncFab in amrex/Src/Base/AMReX_PhysBCFunct.H
// We can't get around this so instead we create an incflo object
//    and use that to access the quantities that aren't passed here.
inline void VelFillBox(Box const& bx, Array4<amrex::Real> const& dest, 
                       const int dcomp, const int numcomp,
                       GeometryData const& geom, const Real time_in, 
                       const BCRec* bcr, 
                       const int bcomp, const int orig_comp)
{
    if (dcomp != 0)
         amrex::Abort("Must have dcomp = 0 in VelFillBox");
    
    if (numcomp != 3)
         amrex::Abort("Must have numcomp = 3 in VelFillBox");
    

    const Box& domain = geom.Domain();

    // This is a bit hack-y but does get us the right level 
    int lev = 0;
    while(lev < 20)
    {
       const Geometry& lev_geom = incflo_for_fillpatching->get_geom_ref(lev);
       if (domain.length()[0] == (lev_geom.Domain()).length()[0]) 
       {
         break;
       }
       lev++;
    }
    if (lev == 20)
        amrex::Abort("Reached lev = 20 in VelFillBox...");

    // We are hard-wiring this fillpatch routine to define the Dirichlet values
    //    at the faces (not the ghost cell center)
    int extrap_dir_bcs = 1;

    // We only do this to make it not const
    Real time = time_in;

    const int* bc_ilo_ptr = incflo_for_fillpatching->get_bc_ilo_ptr(lev);
    const int* bc_ihi_ptr = incflo_for_fillpatching->get_bc_ihi_ptr(lev);
    const int* bc_jlo_ptr = incflo_for_fillpatching->get_bc_jlo_ptr(lev);
    const int* bc_jhi_ptr = incflo_for_fillpatching->get_bc_jhi_ptr(lev);
    const int* bc_klo_ptr = incflo_for_fillpatching->get_bc_klo_ptr(lev);
    const int* bc_khi_ptr = incflo_for_fillpatching->get_bc_khi_ptr(lev);

    int nghost = incflo_for_fillpatching->get_nghost();
    int probtype = incflo_for_fillpatching->get_probtype();

    FArrayBox dest_fab(dest);

    set_velocity_bcs(&time, 
                     dest_fab.dataPtr(), dest_fab.loVect(), dest_fab.hiVect(),
                     bc_ilo_ptr, bc_ihi_ptr, 
                     bc_jlo_ptr, bc_jhi_ptr, 
                     bc_klo_ptr, bc_khi_ptr, 
                     domain.loVect(), domain.hiVect(),
                     &nghost, &extrap_dir_bcs, &probtype);
}

// Compute a new multifab by copying array from valid region and filling ghost cells
// works for single level and 2-level cases (fill fine grid ghost by interpolating from coarse)
void
incflo::FillPatchVel(int lev, Real time, MultiFab& mf, int icomp, int ncomp)
{
    // There aren't used for anything but need to be defined for the function call
    Vector<BCRec> bcs(3);

    // Hack so that ghost cells are not undefined
    mf.setDomainBndry(boundary_val, geom[lev]);

    if (lev == 0)
    {
        Vector<MultiFab*> smf;
        Vector<Real> stime;
        GetDataVel(0, time, smf, stime);

        CpuBndryFuncFab bfunc(VelFillBox);
        PhysBCFunct<CpuBndryFuncFab> physbc(geom[lev], bcs, bfunc);
        amrex::FillPatchSingleLevel(mf, time, smf, stime, 0, icomp, ncomp,
                                    geom[lev], physbc, 0);
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
                                  0, icomp, ncomp, geom[lev-1], geom[lev],
                                  cphysbc, 0, fphysbc, 0,
                                  refRatio(lev-1), mapper, bcs, 0);

    }
}

// utility to copy in data from phi_old and/or phi_new into another multifab
void
incflo::GetDataVel(int lev, Real time, Vector<MultiFab*>& data, Vector<Real>& datatime)
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

void
incflo::FillPatchDensity(int lev, Real time, MultiFab& mf, int icomp, int ncomp)
{
    // There aren't used for anything but need to be defined for the function call
    Vector<BCRec> bcs(3);

    // Hack so that ghost cells are not undefined
    mf.setDomainBndry(boundary_val, geom[lev]);

    if (lev == 0)
    {
        Vector<MultiFab*> smf;
        Vector<Real> stime;
        GetDataDensity(0, time, smf, stime);

        CpuBndryFuncFab bfunc(DensityFillBox);
        PhysBCFunct<CpuBndryFuncFab> physbc(geom[lev], bcs, bfunc);
        amrex::FillPatchSingleLevel(mf, time, smf, stime, 0, icomp, ncomp,
                                    geom[lev], physbc, 0);
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
                                  0, icomp, ncomp, geom[lev-1], geom[lev],
                                  cphysbc, 0, fphysbc, 0,
                                  refRatio(lev-1), mapper, bcs, 0);

    }
}

// utility to copy in data from phi_old and/or phi_new into another multifab
void
incflo::GetDataDensity(int lev, Real time, Vector<MultiFab*>& data, Vector<Real>& datatime)
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
        data.push_back(density[lev].get());
        datatime.push_back(t_old[lev]);
        datatime.push_back(t_new[lev]);
    }
}

void
incflo::FillPatchTracer(int lev, Real time, MultiFab& mf, int icomp, int ncomp)
{
    // There aren't used for anything but need to be defined for the function call
    Vector<BCRec> bcs(3);

    // Hack so that ghost cells are not undefined
    mf.setDomainBndry(boundary_val, geom[lev]);

    if (lev == 0)
    {
        Vector<MultiFab*> smf;
        Vector<Real> stime;
        GetDataTracer(0, time, smf, stime);

        CpuBndryFuncFab bfunc(TracerFillBox);
        PhysBCFunct<CpuBndryFuncFab> physbc(geom[lev], bcs, bfunc);
        amrex::FillPatchSingleLevel(mf, time, smf, stime, 0, icomp, ncomp,
                                    geom[lev], physbc, 0);
    }
    else
    {
        Vector<MultiFab*> cmf, fmf;
        Vector<Real> ctime, ftime;
        GetDataTracer(lev-1, time, cmf, ctime);
        GetDataTracer(lev  , time, fmf, ftime);

        CpuBndryFuncFab bfunc(TracerFillBox);
        PhysBCFunct<CpuBndryFuncFab> cphysbc(geom[lev-1],bcs,bfunc);
        PhysBCFunct<CpuBndryFuncFab> fphysbc(geom[lev  ],bcs,bfunc);

        Interpolater* mapper = &cell_cons_interp;

        amrex::FillPatchTwoLevels(mf, time, cmf, ctime, fmf, ftime,
                                  0, icomp, ncomp, geom[lev-1], geom[lev],
                                  cphysbc, 0, fphysbc, 0,
                                  refRatio(lev-1), mapper, bcs, 0);

    }
}

// utility to copy in data from phi_old and/or phi_new into another multifab
void
incflo::GetDataTracer(int lev, Real time, Vector<MultiFab*>& data, Vector<Real>& datatime)
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
        data.push_back(tracer[lev].get());
        datatime.push_back(t_old[lev]);
        datatime.push_back(t_new[lev]);
    }
}

