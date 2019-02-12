#include <AMReX_EB2.H>
#include <AMReX_EB2_IF_Plane.H>
#include <AMReX_EB2_IF_Union.H>

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
    // Check whether there are walls present
    bool has_walls = false;
    std::unique_ptr<UnionListIF<EB2::PlaneIF>> impfunc_walls = get_real_walls(has_walls);

    // Settings for building index space
    int max_level_here = 0;
	int max_coarsening_level = 100;
    EBSupport m_eb_support_level = EBSupport::full;

    // Generate GeometryShop and build index space
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

    // Make the EBFabFactory
    for(int lev = 0; lev <= max_level; lev++)
    {
        const EB2::Level& eb_is_lev = eb_is.getLevel(geom[lev]);
        eb_level = &eb_is_lev;
        ebfactory[lev].reset(new EBFArrayBoxFactory(*eb_level, 
                                                    geom[lev], 
                                                    grids[lev], 
                                                    dmap[lev],
                                                    {m_eb_basic_grow_cells, 
                                                    m_eb_volume_grow_cells, 
                                                    m_eb_full_grow_cells},
                                                    m_eb_support_level));
    }
}
