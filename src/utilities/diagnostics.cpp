#include <incflo.H>

using namespace amrex;

//
// Print maximum values (useful for tracking evolution)
//
void incflo::PrintMaxValues(Real time_in)
{
    return; // xxxxx TODO

#if 0
    ComputeDivU(time_in);

    for(int lev = 0; lev <= finest_level; lev++)
    {
        amrex::Print() << "Level " << lev << std::endl;
        PrintMaxVel(lev);
        PrintMaxGp(lev);
    }
    amrex::Print() << std::endl;
#endif
}

//
// Print the maximum values of the velocity components and velocity divergence
//
void incflo::PrintMaxVel(int lev)
{
#if 0
    // xxxxx TODO
    amrex::Print() << "max(abs(u/v/w/divu))  = "
                   << Norm(vel , lev, 0, 0) << "  "
		   << Norm(vel , lev, 1, 0) << "  "
                   << Norm(vel , lev, 2, 0) << "  "
                   << Norm(divu, lev, 0, 0) << "  " << std::endl;
    if (ntrac > 0)
    {
       for (int i = 0; i < ntrac; i++)
          amrex::Print() << "max tracer" << i << " = " << Norm(tracer,lev,i,0) << std::endl;;
    }
#endif
}

//
// Print the maximum values of the pressure gradient components and pressure
//
void incflo::PrintMaxGp(int lev)
{
#if 0
    // xxxxx TODO
    amrex::Print() << "max(abs(gpx/gpy/gpz/p))  = "
                   << Norm(gp, lev, 0, 0) << "  "
		   << Norm(gp, lev, 1, 0) << "  "
                   << Norm(gp, lev, 2, 0) << "  "
		   << Norm(p , lev, 0, 0) << "  " << std::endl;
#endif
}

void incflo::CheckForNans(int lev)
{
#if 0
    // xxxxx
    bool ro_has_nans = density[lev]->contains_nan(0);
    bool ug_has_nans = vel[lev]->contains_nan(0);
    bool vg_has_nans = vel[lev]->contains_nan(1);
    bool wg_has_nans = vel[lev]->contains_nan(2);
    bool pg_has_nans = p[lev]->contains_nan(0);

    if (ro_has_nans)
	amrex::Print() << "WARNING: ro contains NaNs!!!";

    if (ug_has_nans)
	amrex::Print() << "WARNING: u contains NaNs!!!";

    if(vg_has_nans)
	amrex::Print() << "WARNING: v contains NaNs!!!";

    if(wg_has_nans)
	amrex::Print() << "WARNING: w contains NaNs!!!";

    if(pg_has_nans)
	amrex::Print() << "WARNING: p contains NaNs!!!";
#endif
}

Real
incflo::volWgtSum (int lev, const MultiFab& mf, int comp, bool local)
{
    return 0;
#if 0
    // xxxxx
    BL_PROFILE("incflo::volWgtSum()");

#ifdef AMREX_USE_EB
    const MultiFab* volfrac =  &(ebfactory[lev]->getVolFrac());
#endif

#ifdef AMREX_USE_CUDA
    bool switch_GPU_launch(false);

    if(Gpu::notInLaunchRegion())
    {
      switch_GPU_launch = true;
      Gpu::setLaunchRegion(true);
    }
#endif

#ifdef AMREX_USE_EB
    Real sum = amrex::ReduceSum(mf, *volfrac, 0,
        [comp] AMREX_GPU_HOST_DEVICE (Box const & bx,
                                      FArrayBox const & rho_fab,
                                      FArrayBox const & volfrac_fab)
        {
          Real dm = 0.0;
          const auto rho = rho_fab.const_array();
          const auto vfrc = volfrac_fab.const_array();

          amrex::Loop(bx, [rho,vfrc,comp,&dm] (int i, int j, int k) noexcept
              { dm += rho(i,j,k,comp) * vfrc(i,j,k); });

          return dm;
        });
#else
    Real sum = amrex::ReduceSum(mf, 0,
        [comp] AMREX_GPU_HOST_DEVICE (Box const & bx,
                                      FArrayBox const & rho_fab)
        {
          Real dm = 0.0;
          const auto rho = rho_fab.const_array();

          amrex::Loop(bx, [rho,comp,&dm] (int i, int j, int k) noexcept
              { dm += rho(i,j,k,comp); });

          return dm;
        });
#endif

#ifdef AMREX_USE_CUDA
    if(switch_GPU_launch)
      Gpu::setLaunchRegion(false);
#endif

    if (!local)
        ParallelDescriptor::ReduceRealSum(sum);

    return sum;
#endif
}

