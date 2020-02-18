#include <incflo.H>

using namespace amrex;

//
// Print maximum values (useful for tracking evolution)
//
void incflo::PrintMaxValues(Real time_in)
{
    
#if AMREX_USE_EB
    amrex::Print << "Warning PrintMaxValues() does not work with EB yet" << std::endl;
    return; // xxxxx TODO
#endif
    
//    ComputeDivU(time_in);

    for(int lev = 0; lev <= finest_level; lev++)
    {
        amrex::Print() << "Level " << lev << std::endl;
        PrintMaxVel(lev);
        PrintMaxGp(lev);
    }
    amrex::Print() << std::endl;

}

//
// Print the maximum values of the velocity components and velocity divergence
//
void incflo::PrintMaxVel(int lev)
{

    LevelData &ld = *m_leveldata[lev];
    
    amrex::Print() << "max(abs(u/v/w))  = "
                   << ld.velocity.norm0(0) << "  "
                   << ld.velocity.norm0(1)  << "  "
                   << ld.velocity.norm0(2)  << "  "
                   << std::endl;
    if (m_ntrac > 0)
    {
       for (int i = 0; i < m_ntrac; i++)
          amrex::Print() << "max tracer" << i << " = " << ld.tracer.norm0(i) << std::endl;;
    }

}

//
// Print the maximum values of the pressure gradient components and pressure
//
void incflo::PrintMaxGp(int lev)
{

    LevelData &ld = *m_leveldata[lev];

    amrex::Print() << "max(abs(gpx/gpy/gpz/p))  = "
                   << ld.gp.norm0(0) << "  "
                   << ld.gp.norm0(1) << "  "
                   << ld.gp.norm0(2) << "  "
		           << ld.p.norm0(0) << "  " << std::endl;

}

void incflo::CheckForNans(int lev)
{

    LevelData &ld = *m_leveldata[lev];

    bool ro_has_nans = ld.density.contains_nan(0);
    bool ug_has_nans = ld.velocity.contains_nan(0);
    bool vg_has_nans = ld.velocity.contains_nan(1);
    bool wg_has_nans = ld.velocity.contains_nan(2);
    bool pg_has_nans = ld.p.contains_nan(0);

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

