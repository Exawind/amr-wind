#include <AMReX_EB2.H>
#include <AMReX_EB2_IF_Sphere.H>
#include <AMReX_ParmParse.H>

#include <algorithm>
#include <incflo.H>

/********************************************************************************
 *                                                                              *
 * Function to create a simple sphere EB.                                     *
 *                                                                              *
 ********************************************************************************/
void incflo::make_eb_sphere()
{
    ParmParse pp("sphere");

    int max_level_here = 0;

    /****************************************************************************
     * Get sphere information from inputs file.                               *
     ***************************************************************************/
	bool inside = true;

    Real radius = 0.0002;

    int direction = 0;
    Vector<Real> centervec(3);

    pp.query("internal_flow", inside);

    pp.query("radius", radius);
    pp.getarr("center", centervec, 0, 3);
    Array<Real, 3> center = {centervec[0], centervec[1], centervec[2]};

	amrex::Print() << " " << std::endl;
	amrex::Print() << " Internal Flow: " << inside << std::endl;
	amrex::Print() << " # Ghost Cells: " << nghost << std::endl;
	amrex::Print() << " Radius:    " << radius << std::endl;
	amrex::Print() << " Center:    " << center[0] << ", " << center[1] << ", " << center[2]
				   << std::endl;

    /****************************************************************************
     *                                                                          *
     * Build standard EB Factories                                              *
     *                                                                          *
     ***************************************************************************/

    // set up ebfactory
    EBSupport m_eb_support_level = EBSupport::full;

    // Create the sphere
    amrex::Print() << "Building the sphere geometry ..." << std::endl;

    EB2::SphereIF my_sphere(radius, center, inside);

    auto gshop = EB2::makeShop(my_sphere);
	int max_coarsening_level = 100;
	EB2::Build(gshop, geom.back(), max_level_here, max_level_here + max_coarsening_level);
    const EB2::IndexSpace& ebis = EB2::IndexSpace::top();

    amrex::Print() << "Done building the sphere geometry" << std::endl;

    /****************************************************************************
    *                                                                           *
    * THIS FILLS FLUID EBFACTORY                                                *
    *                                                                           *
    *****************************************************************************/
    amrex::Print() << "Now  making the fluid ebfactory ..." << std::endl;

    for(int lev = 0; lev <= max_level; lev++)
    {
        const EB2::Level& ebis_lev = ebis.getLevel(geom[lev]);
        eb_level_fluid = &ebis_lev;
        ebfactory[lev].reset(new EBFArrayBoxFactory(
            *eb_level_fluid,
            geom[lev],
            grids[lev],
            dmap[lev],
            {m_eb_basic_grow_cells, m_eb_volume_grow_cells, m_eb_full_grow_cells},
            m_eb_support_level));
    }

    amrex::Print() << "Done making the fluid ebfactory ..." << std::endl;
}
