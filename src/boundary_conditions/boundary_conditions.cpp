#include <AMReX_Array.H>
#include <AMReX_BC_TYPES.H>
#include <AMReX_BLassert.H>
#include <AMReX_Box.H>
#include <AMReX_EBMultiFabUtil.H>
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


//
// Fill the BCs for velocity only
//
void incflo::FillVelocityBC(Real time, int extrap_dir_bcs)
{
    BL_PROFILE("incflo::FillVelocityBC()");

    for(int lev = 0; lev <= finest_level; lev++)
    {
        Box domain(geom[lev].Domain());

        // Hack so that ghost cells are not undefined
        vel[lev]->setDomainBndry(boundary_val, geom[lev]);

        vel[lev]->FillBoundary(geom[lev].periodicity());
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for(MFIter mfi(*vel[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            set_velocity_bcs(&time, 
                             BL_TO_FORTRAN_ANYD((*vel[lev])[mfi]),
                             bc_ilo[lev]->dataPtr(), bc_ihi[lev]->dataPtr(),
                             bc_jlo[lev]->dataPtr(), bc_jhi[lev]->dataPtr(),
                             bc_klo[lev]->dataPtr(), bc_khi[lev]->dataPtr(),
                             domain.loVect(), domain.hiVect(),
                             &nghost, &extrap_dir_bcs, &probtype);
        }
        EB_set_covered(*vel[lev], covered_val);
        
        // Do this after as well as before to pick up terms that got updated in the call above
        vel[lev]->FillBoundary(geom[lev].periodicity());
    }
}

void incflo::FillScalarBC()
{
    BL_PROFILE("incflo:FillScalarBC()");

    for(int lev = 0; lev <= finest_level; lev++)
    {
        Box domain(geom[lev].Domain());
        
        // Hack so that ghost cells are not undefined
         ro[lev]->setDomainBndry(boundary_val, geom[lev]);
        eta[lev]->setDomainBndry(boundary_val, geom[lev]);

        // Impose periodic BCs at domain boundaries and fine-fine copies in the interior
         ro[lev]->FillBoundary(geom[lev].periodicity());
        eta[lev]->FillBoundary(geom[lev].periodicity());

        // Fill all cell-centered arrays with first-order extrapolation at domain boundaries
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for(MFIter mfi(*ro[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            // Density
            fill_bc0(BL_TO_FORTRAN_ANYD((*ro[lev])[mfi]),
                     bc_ilo[lev]->dataPtr(), bc_ihi[lev]->dataPtr(),
                     bc_jlo[lev]->dataPtr(), bc_jhi[lev]->dataPtr(),
                     bc_klo[lev]->dataPtr(), bc_khi[lev]->dataPtr(),
                     domain.loVect(), domain.hiVect(),
                     &nghost);

            // Viscosity
            fill_bc0(BL_TO_FORTRAN_ANYD((*eta[lev])[mfi]),
                     bc_ilo[lev]->dataPtr(), bc_ihi[lev]->dataPtr(),
                     bc_jlo[lev]->dataPtr(), bc_jhi[lev]->dataPtr(),
                     bc_klo[lev]->dataPtr(), bc_khi[lev]->dataPtr(),
                     domain.loVect(), domain.hiVect(),
                     &nghost);
        }
    }
}

void incflo::GetInputBCs()
{
    // Extracts all walls from the inputs file
    int cyclic;

    cyclic = geom[0].isPeriodic(0) ? 1 : 0;
    SetInputBCs("xlo", 1, cyclic, geom[0].ProbLo(0));
    SetInputBCs("xhi", 2, cyclic, geom[0].ProbHi(0));

    cyclic = geom[0].isPeriodic(1) ? 1 : 0;
    SetInputBCs("ylo", 3, cyclic, geom[0].ProbLo(1));
    SetInputBCs("yhi", 4, cyclic, geom[0].ProbHi(1));

    cyclic = geom[0].isPeriodic(2) ? 1 : 0;
    SetInputBCs("zlo", 5, cyclic, geom[0].ProbLo(2));
    SetInputBCs("zhi", 6, cyclic, geom[0].ProbHi(2));
}

void incflo::SetInputBCs(const std::string bcID, const int index,
                           const int cyclic, const Real domloc) 
{
    const int und_  =   0;
    const int pinf_ =  10;
    const int pout_ =  11;
    const int minf_ =  20;
    const int nsw_  = 100;

    // Default a BC to undefined.
    int itype = und_;

    int direction = 0;
    Real pressure = -1.0;
    Vector<Real> velocity(3, 0.0);
    Real location = domloc;

    std::string bc_type = "null";

    ParmParse pp(bcID);

    pp.query("type", bc_type);

    if(bc_type == "pressure_inflow"  || bc_type == "pi" ||
              bc_type == "PRESSURE_INFLOW"  || bc_type == "PI" ) {

      amrex::Print() << bcID <<" set to pressure inflow. "  << std::endl;
      itype = pinf_;

      pp.get("pressure", pressure);

    } else if(bc_type == "pressure_outflow" || bc_type == "po" ||
              bc_type == "PRESSURE_OUTFLOW" || bc_type == "PO" ) {

      amrex::Print() << bcID <<" set to pressure outflow. "  << std::endl;
      itype = pout_;

      pp.get("pressure", pressure);


    } else if (bc_type == "mass_inflow"     || bc_type == "mi" ||
               bc_type == "MASS_INFLOW"     || bc_type == "MI" ) {

      // Flag that this is a mass inflow.
      amrex::Print() << bcID <<" set to mass inflow. "  << std::endl;
      itype = minf_;

      pp.query("pressure", pressure);
      pp.getarr("velocity", velocity, 0, 3);


    } else if (bc_type == "no_slip_wall"    || bc_type == "nsw" ||
               bc_type == "NO_SLIP_WALL"    || bc_type == "NSW" ) {

      // Flag that this is a no-slip wall.
      amrex::Print() << bcID <<" set to no-slip wall. "  << std::endl;
      itype = nsw_;

      pp.queryarr("velocity", velocity, 0, 3);
      pp.query("direction", direction);
      pp.query("location", location);

    }

    if ( cyclic == 1 && itype != und_){
      amrex::Abort("Cannot mix periodic BCs and Wall/Flow BCs.\n");
    }

    const Real* plo = geom[0].ProbLo();
    const Real* phi = geom[0].ProbHi();

    set_bc_mod(&index, &itype, plo, phi,
               &location, &pressure, &velocity[0]);

}
