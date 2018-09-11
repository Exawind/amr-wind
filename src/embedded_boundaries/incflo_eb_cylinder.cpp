#include <AMReX_EB2.H>
#include <AMReX_EB2_IF_Cylinder.H>
#include <AMReX_EB2_IF_Plane.H>
#include <AMReX_EB2_IF_Union.H>

#include <AMReX_ParmParse.H>

//#include <AMReX_VisMF.H>  // amrex::VisMF::Write(MultiFab)
//#include <sstream>

#include <AMReX_EB_levelset.H>
#include <algorithm>
#include <incflo_eb_F.H>
#include <incflo_level.H>

using namespace std;

/********************************************************************************
 *                                                                              *
 * Function to create a simple cylinder EB.                                     *
 *                                                                              *
 * Comments: The cylinder can either be closed or open at the bottom. At this   *
 * time, there is no option to put a "lid" on the top of the cylinder.          *
 *                                                                              *
 * The bottom of the cylinder is only applied to particles. This allows for the *
 * fluid to not have an EB surface generated that would block gas flow from a   *
 * mass inflow boundary. A more general implementation would check the incflo.dat *
 * file for inflows and outflows and correctly cap the top and bottom of the    *
 * cylinder.                                                                    *
 *                                                                              *
 ********************************************************************************/
void incflo_level::make_eb_cylinder(int lev)
{
	ParmParse pp("cylinder");

	int max_level_here = 0;

	/****************************************************************************
     * Get cylinder information from inputs file.                               *
     ****************************************************************************/
	bool inside = true;
	bool close_bottom = true;
	Real offset = 1.0e-8;

	Real radius = 0.0002;
	Real height = 0.0080;

	int direction = 0;
	Vector<Real> centervec(3);

	pp.query("internal_flow", inside);
	pp.query("closed_bottom", close_bottom);
	pp.query("bottom_offset", offset);

	pp.query("radius", radius);
	pp.query("height", height);
	pp.query("direction", direction);
	pp.getarr("center", centervec, 0, 3);
	Array<Real, 3> center = {centervec[0], centervec[1], centervec[2]};

	/****************************************************************************
     *                                                                          *
     * Build standard EB Factories                                              *
     *                                                                          *
     ****************************************************************************/

	// set up ebfactory
	int m_eb_basic_grow_cells = nghost;
	int m_eb_volume_grow_cells = nghost;
	int m_eb_full_grow_cells = nghost;
	EBSupport m_eb_support_level = EBSupport::full;

	amrex::Print() << " " << std::endl;
	amrex::Print() << " Internal Flow: " << inside << std::endl;
	amrex::Print() << " Radius:    " << radius << std::endl;
	amrex::Print() << " Height:    " << height << std::endl;
	amrex::Print() << " Offset:    " << offset << std::endl;
	amrex::Print() << " Direction: " << direction << std::endl;
	amrex::Print() << " Center:    " << center[0] << ", " << center[1] << ", " << center[2]
				   << std::endl;

	// Create the cylinder
	amrex::Print() << "Building the cylinder (side wall) geometry ..." << std::endl;

	EB2::CylinderIF my_cyl(radius, height, direction, center, inside);

	auto gshop_cyl = EB2::makeShop(my_cyl);
	int max_coarsening_level = 100;
	EB2::Build(gshop_cyl, geom.back(), max_level_here, max_level_here + max_coarsening_level);

	const EB2::IndexSpace& ebis_cyl = EB2::IndexSpace::top();
	const EB2::Level& ebis_lev_cyl = ebis_cyl.getLevel(geom[lev]);

	amrex::Print() << "Done building the cylinder geometry" << std::endl;

	/****************************************************************************
    *                                                                           *
    * THIS FILLS FLUID EBFACTORY                                                *
    *                                                                           *
    *****************************************************************************/

	amrex::Print() << "Now  making the fluid ebfactory ..." << std::endl;

	eb_level_fluid = &ebis_lev_cyl;

	ebfactory[lev].reset(new EBFArrayBoxFactory(
		*eb_level_fluid,
		geom[lev],
		grids[lev],
		dmap[lev],
		{m_eb_basic_grow_cells, m_eb_volume_grow_cells, m_eb_full_grow_cells},
		m_eb_support_level));

	amrex::Print() << "Done making the fluid ebfactory ..." << std::endl;
}
