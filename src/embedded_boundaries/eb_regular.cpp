#include <AMReX_EB2.H>
#include <AMReX_EB2_IF_Plane.H>
#include <AMReX_EB2_IF_Union.H>
#include <AMReX_EB_levelset.H>

#include <eb_if.H>
#include <embedded_boundaries_F.H>
#include <incflo.H>

/********************************************************************************
 *                                                                              *
 * Placeholder: create a simulation box _without_ EB walls.                     *
 *                                                                              *
 ********************************************************************************/
void incflo::make_eb_regular()
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

	std::unique_ptr<MultiFab> mf_impfunc;

    bool has_walls = false;
    std::unique_ptr<UnionListIF<EB2::PlaneIF>> impfunc_walls = get_real_walls(has_walls);

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

    for(int lev = 0; lev <= max_level; lev++)
    {
        eb_level_fluid = &eb_is.getLevel(geom[lev]);
        ebfactory[lev].reset(new EBFArrayBoxFactory(*eb_level_fluid,
                                                    geom[lev],
                                                    grids[lev],
                                                    dmap[lev],
                                                    {
                                                        m_eb_basic_grow_cells, 
                                                        m_eb_volume_grow_cells, 
                                                        m_eb_full_grow_cells
                                                    },
                                                    m_eb_support_level));
    }
}
