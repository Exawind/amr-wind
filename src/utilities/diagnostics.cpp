#include <incflo.H>

using namespace amrex;

//
// Print maximum values (useful for tracking evolution)
//
void incflo::PrintMaxValues(Real /* time_in */)
{
    
    
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
