#include <AMReX_EB2.H>
#include <AMReX_EB2_IF_Plane.H>
#include <AMReX_EB2_IF_Union.H>

#include <incflo_eb_if.H>

#include <AMReX_EB_levelset.H>

#include <incflo_eb_F.H>
#include <incflo_level.H>

/********************************************************************************
 *                                                                              *
 * Placeholder: create a simulation box _without_ EB walls.                     *
 *                                                                              *
 ********************************************************************************/
void incflo_level::make_eb_regular(int lev)
{
	int max_level_here = 0;

	/****************************************************************************
     *                                                                          *
     * Build standard EB Factories                                              *
     *                                                                          *
     ****************************************************************************/

	// set up ebfactory

	EBSupport m_eb_support_level = EBSupport::full;

	int max_coarsening_level = 100;

	amrex::Print() << " " << std::endl;
	amrex::Print() << "Now making the ebfactories ..." << std::endl;

	// If filling level-set: this is used to store the implicit function (due to
	// any walls defined in incflo.dat). It is filled while after EB2::Build.
	// NOTE: this pointer is undefined if and of:
	//     * ! solve_dem
	//     * levelset__restart
	//     * ! has_walls
	// are true
	std::unique_ptr<MultiFab> mf_impfunc;

	bool has_walls = false;
	std::unique_ptr<UnionListIF<EB2::PlaneIF>> impfunc_walls = get_real_walls(lev, has_walls);

	if(has_walls)
	{
		auto gshop = EB2::makeShop(*impfunc_walls);
		EB2::Build(gshop, geom.back(), max_level_here, max_level_here + max_coarsening_level);
	}
	else
	{
		EB2::AllRegularIF my_regular;
		auto gshop = EB2::makeShop(my_regular);
		EB2::Build(gshop, geom.back(), max_level_here, max_level_here + max_coarsening_level);
	}

	const EB2::IndexSpace& eb_is = EB2::IndexSpace::top();
	eb_level_fluid = &eb_is.getLevel(geom[lev]);

	// Do _not_ fill level-set with AllRegularIF => if there are no walls, then
	// the level-set function is just huge(amrex_real) => this flag is set to
	// true iff there are walls.
	has_walls = false;

	ebfactory[lev].reset(new EBFArrayBoxFactory(
		*eb_level_fluid,
		geom[lev],
		grids[lev],
		dmap[lev],
		{m_eb_basic_grow_cells, m_eb_volume_grow_cells, m_eb_full_grow_cells},
		m_eb_support_level));
}
