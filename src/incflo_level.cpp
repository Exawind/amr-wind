#include <AMReX_ParmParse.H>

#include <AMReX_BC_TYPES.H>
#include <AMReX_Box.H>
#include <AMReX_EBMultiFabUtil.H>

#include <incflo_level.H>

// Initiate vars which cannot be initiated in header
Vector<Real> incflo_level::gravity(3, 0.);
std::string incflo_level::load_balance_type = "FixedSize";
std::string incflo_level::knapsack_weight_type = "RunTimeCosts";

// Define unit vectors for easily convert indeces
amrex::IntVect incflo_level::e_x(1, 0, 0);
amrex::IntVect incflo_level::e_y(0, 1, 0);
amrex::IntVect incflo_level::e_z(0, 0, 1);

EBSupport incflo_level::m_eb_support_level = EBSupport::full;

incflo_level::~incflo_level(){};

incflo_level::incflo_level()
{
// Geometry on all levels has just been defined in the AmrCore constructor

// No valid BoxArray and DistributionMapping have been defined.
// But the arrays for them have been resized.

#if 0
    int nlevs_max = maxLevel() + 1;
    istep.resize(nlevs_max, 0);
    nsubsteps.resize(nlevs_max, 1);
    for (int lev = 1; lev <= maxLevel(); ++lev) 
        nsubsteps[lev] = MaxRefRatio(lev-1);
#endif
}

//
// Subroutine to compute norm0 of EB multifab
//
Real incflo_level::incflo_norm0(const Vector<std::unique_ptr<MultiFab>>& mf, 
                                int lev, int comp)
{
	MultiFab mf_tmp(mf[lev]->boxArray(),
					mf[lev]->DistributionMap(),
                    mf[lev]->nComp(), 
                    0, MFInfo(), *ebfactory[lev]);

	MultiFab::Copy(mf_tmp, *mf[lev], comp, comp, 1, 0);
	EB_set_covered(mf_tmp, 0.0);

	return mf_tmp.norm0(comp);
}

Real incflo_level::incflo_norm0(MultiFab& mf, int lev, int comp)
{
    MultiFab mf_tmp(mf.boxArray(), 
                    mf.DistributionMap(), 
                    mf.nComp(), 
                    0, MFInfo(), *ebfactory[lev]);

	MultiFab::Copy(mf_tmp, mf, comp, comp, 1, 0);
	EB_set_covered(mf_tmp, 0.0);

	return mf_tmp.norm0(comp);
}

//
// Subroutine to compute norm1 of EB multifab
//
Real incflo_level::incflo_norm1(const Vector<std::unique_ptr<MultiFab>>& mf, 
                                int lev, int comp)
{
	MultiFab mf_tmp(mf[lev]->boxArray(),
					mf[lev]->DistributionMap(),
					mf[lev]->nComp(),
                    0, MFInfo(), *ebfactory[lev]);

	MultiFab::Copy(mf_tmp, *mf[lev], comp, comp, 1, 0);
	EB_set_covered(mf_tmp, 0.0);

	return mf_tmp.norm1(comp, geom[lev].periodicity());
}

Real incflo_level::incflo_norm1(MultiFab& mf, int lev, int comp)
{
	MultiFab mf_tmp(mf.boxArray(), 
                    mf.DistributionMap(), 
                    mf.nComp(), 
                    0, MFInfo(), *ebfactory[lev]);

	MultiFab::Copy(mf_tmp, mf, comp, comp, 1, 0);
	EB_set_covered(mf_tmp, 0.0);

	return mf_tmp.norm1(comp, geom[lev].periodicity());
}

//
// Print the maximum values of the velocity components
//
void incflo_level::incflo_print_max_vel(int lev)
{
	amrex::Print() << "max(abs(u/v/w/p))  = " 
                   << incflo_norm0(vel, lev, 0) << "  "
				   << incflo_norm0(vel, lev, 1) << "  " 
                   << incflo_norm0(vel, lev, 2) << "  "
				   << incflo_norm0(p, lev, 0) << "  " << std::endl;
}

void incflo_level::check_for_nans(int lev)
{
	bool ug_has_nans = vel[lev]->contains_nan(0);
	bool vg_has_nans = vel[lev]->contains_nan(1);
	bool wg_has_nans = vel[lev]->contains_nan(2);
	bool pg_has_nans = p[lev]->contains_nan(0);

	if(ug_has_nans)
		amrex::Print() << "WARNING: u contains NaNs!!!";

	if(vg_has_nans)
		amrex::Print() << "WARNING: v contains NaNs!!!";

	if(wg_has_nans)
		amrex::Print() << "WARNING: w contains NaNs!!!";

	if(pg_has_nans)
		amrex::Print() << "WARNING: p contains NaNs!!!";
}


void incflo_level::Regrid(int base_lev)
{
	BL_PROFILE_REGION_START("incflo::Regrid()");

	if(load_balance_type == "KnapSack")
	{
		AMREX_ALWAYS_ASSERT(fluid_cost[0] != nullptr);

		if(ParallelDescriptor::NProcs() == 1)
			return;
		{
			MultiFab costs(grids[base_lev], dmap[base_lev], 1, 0);
			costs.setVal(0.0);
			costs.plus(*fluid_cost[base_lev], 0, 1, 0);

			DistributionMapping newdm = DistributionMapping::makeKnapSack(costs);

			bool dm_changed = (newdm != dmap[base_lev]);

			SetDistributionMap(base_lev, newdm);

			if(dm_changed)
				RegridArrays(base_lev);

			fluid_cost[base_lev].reset(new MultiFab(grids[base_lev], newdm, 1, 0));
			fluid_cost[base_lev]->setVal(0.0);

			incflo_set_bc0(base_lev);

			if(ebfactory[base_lev])
			{
				ebfactory[base_lev].reset(new EBFArrayBoxFactory(
					*eb_level_fluid,
					geom[base_lev],
					grids[base_lev],
					dmap[base_lev],
					{m_eb_basic_grow_cells, m_eb_volume_grow_cells, m_eb_full_grow_cells},
					m_eb_support_level));
			}
		}
	}
	BL_PROFILE_REGION_STOP("incflo::Regrid()");
}
