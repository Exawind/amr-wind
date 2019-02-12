#include <AMReX_EB2.H>
#include <AMReX_EB2_IF_Cylinder.H>
#include <AMReX_ParmParse.H>

#include <algorithm>
#include <incflo.H>

/********************************************************************************
 *                                                                              *
 * Function to create a simple cylinder EB.                                     *
 *                                                                              *
 ********************************************************************************/
void incflo::make_eb_cylinder()
{
    // Initialise cylinder parameters
	bool inside = true;
    Real radius = 0.0002;
    int direction = 0;
    Vector<Real> centervec(3);

    // Get cylinder information from inputs file.                               *
    ParmParse pp("cylinder");

    pp.query("internal_flow", inside);
    pp.query("radius", radius);
    pp.query("direction", direction);
    pp.getarr("center", centervec, 0, 3);
    Array<Real, 3> center = {centervec[0], centervec[1], centervec[2]};

    // Print info about cylinder
	amrex::Print() << " " << std::endl;
	amrex::Print() << " Internal Flow: " << inside << std::endl;
	amrex::Print() << " Radius:    " << radius << std::endl;
	amrex::Print() << " Direction: " << direction << std::endl;
	amrex::Print() << " Center:    " << center[0] << ", " << center[1] << ", " << center[2]
				   << std::endl;

    // Build the Cylinder implficit function representing the curved walls     
    EB2::CylinderIF my_cyl(radius, direction, center, inside);

    // Generate GeometryShop
    auto gshop = EB2::makeShop(my_cyl);

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
