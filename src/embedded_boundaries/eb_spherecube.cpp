#include <AMReX_EB2.H>
#include <AMReX_EB2_IF.H>
#include <AMReX_ParmParse.H>

#include <algorithm>
#include <incflo.H>

/********************************************************************************
 *                                                                              *
 * Function to create a simple sphereCUBE EB.                                     *
 *                                                                              *
 ********************************************************************************/
void incflo::make_eb_spherecube()
{
	amrex::Print() << " " << std::endl;
	amrex::Print() << " SPERECUUUUBE: " << std::endl;

    // Build the sphere implicit function 
    EB2::SphereIF sphere(0.5, {1.8, 1.8, 2.8}, false);
    EB2::BoxIF cube({1.85, 1.85, 2.85}, {2.5, 2.5, 3.5}, false);
    auto cubesphere = EB2::makeUnion(sphere, cube);

    // Generate GeometryShop
    auto gshop = EB2::makeShop(cubesphere);

    // Build index space
    int max_level_here = 0;
	int max_coarsening_level = 100;
    EBSupport m_eb_support_level = EBSupport::full;
	EB2::Build(gshop, geom.back(), max_level_here, max_level_here + max_coarsening_level);
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
