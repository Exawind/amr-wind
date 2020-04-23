#include <incflo.H>

using namespace amrex;

//
// Print maximum values (useful for tracking evolution)
//
void incflo::PrintMaxValues(Real /* time_in */)
{
    
    BL_PROFILE("amr-wind::incflo::PrintMaxValues")

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
    BL_PROFILE("amr-wind::incflo::PrintMaxVel")
    amrex::Print() << "max(abs(u/v/w))  = "
                   << velocity()(lev).norm0(0) << "  "
                   << velocity()(lev).norm0(1)  << "  "
                   << velocity()(lev).norm0(2)  << "  "
                   << std::endl;

    amrex::Print() << "max tracer = "
                   << tracer()(lev).norm0(0) << std::endl;
       
}

//
// Print the maximum values of the pressure gradient components and pressure
//
void incflo::PrintMaxGp(int lev)
{
    BL_PROFILE("amr-wind::incflo::PrintMaxGp")
    amrex::Print() << "max(abs(gpx/gpy/gpz/p))  = "
                   << grad_p()(lev).norm0(0) << "  "
                   << grad_p()(lev).norm0(1) << "  "
                   << grad_p()(lev).norm0(2) << "  "
                   << pressure()(lev).norm0(0) << "  " << std::endl;

}

void incflo::CheckForNans(int lev)
{
    BL_PROFILE("amr-wind::incflo::CheckForNans")
    bool ro_has_nans = density()(lev).contains_nan(0);
    bool ug_has_nans = velocity()(lev).contains_nan(0);
    bool vg_has_nans = velocity()(lev).contains_nan(1);
    bool wg_has_nans = velocity()(lev).contains_nan(2);
    bool pg_has_nans = pressure()(lev).contains_nan(0);

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
