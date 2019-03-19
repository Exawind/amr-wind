#include <AMReX_AmrCore.H>
#include <AMReX_EBMultiFabUtil.H>

#include <incflo.H>

//
// Subroutine to compute norm of EB multifab
//
Real incflo::Norm(const Vector<std::unique_ptr<MultiFab>>& mf, int lev, int comp, int norm_type)
{
    // Make copy of MF so that we can set values in covered cells to zero. 
	MultiFab mf_tmp(mf[lev]->boxArray(),
					mf[lev]->DistributionMap(),
                    mf[lev]->nComp(),
                    0, MFInfo(), *ebfactory[lev]);

	MultiFab::Copy(mf_tmp, *mf[lev], comp, comp, 1, 0);
	EB_set_covered(mf_tmp, 0.0);
    if(norm_type == 0)
    {
        return mf_tmp.norm0(comp);
    }
    else if(norm_type == 1)
    {
        return mf_tmp.norm1(comp, geom[lev].periodicity());
    }
    else
    {
        amrex::Print() << "Warning: called incflo::Norm() with norm_type not in {0,1}" << std::endl; 
        return -1.0;
    }
}

// 
// Print maximum values (useful for tracking evolution)
void incflo::PrintMaxValues(Real time)
{
        ComputeDivU(time);
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
	amrex::Print() << "max(abs(u/v/w/divu))  = "
                   << Norm(vel, lev, 0, 0) << "  "
				   << Norm(vel, lev, 1, 0) << "  "
                   << Norm(vel, lev, 2, 0) << "  " 
                   << Norm(divu, lev, 0, 0) << "  " << std::endl;
}

//
// Print the maximum values of the pressure gradient components and pressure
//
void incflo::PrintMaxGp(int lev)
{
	amrex::Print() << "max(abs(gpx/gpy/gpz/p))  = "
                   << Norm(gp, lev, 0, 0) << "  "
				   << Norm(gp, lev, 1, 0) << "  "
                   << Norm(gp, lev, 2, 0) << "  "
				   << Norm(p, lev, 0, 0) << "  " << std::endl;
}

void incflo::CheckForNans(int lev)
{
	bool ro_has_nans = ro[lev]->contains_nan(0);
	bool ug_has_nans = vel[lev]->contains_nan(0);
	bool vg_has_nans = vel[lev]->contains_nan(1);
	bool wg_has_nans = vel[lev]->contains_nan(2);
	bool pg_has_nans = p[lev]->contains_nan(0);

	if(ro_has_nans)
    {
		amrex::Print() << "WARNING: ro contains NaNs!!!";
    }

	if(ug_has_nans)
    {
		amrex::Print() << "WARNING: u contains NaNs!!!";
    }

	if(vg_has_nans)
    {
		amrex::Print() << "WARNING: v contains NaNs!!!";
    }

	if(wg_has_nans)
    {
		amrex::Print() << "WARNING: w contains NaNs!!!";
    }

	if(pg_has_nans)
    {
		amrex::Print() << "WARNING: p contains NaNs!!!";
    }
}
