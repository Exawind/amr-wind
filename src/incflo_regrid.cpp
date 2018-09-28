#include <AMReX_ParmParse.H>

#include <AMReX_BC_TYPES.H>
#include <AMReX_Box.H>
#include <incflo_eb_F.H>
#include <incflo_level.H>

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
