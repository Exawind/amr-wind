#include <AMReX_EB2.H>
#include <AMReX_EB2_IF_Cylinder.H>
#include <AMReX_EB2_IF_Union.H>

#include <AMReX_ParmParse.H>

#include <algorithm>
#include <embedded_boundaries_F.H>
#include <incflo.H>

/********************************************************************************
 *                                                                              *
 * Function to create a annular cylinder EB.                                     *
 *                                                                              *
 ********************************************************************************/
void incflo::make_eb_twocylinders()
{
    // Initialise parameters
	int direction1 = 2;
	int direction2 = 2;
	Real radius1 = 0.5;
	Real radius2 = 0.5;
	Vector<Real> centervec1(3);
	Vector<Real> centervec2(3);

    // Get information from inputs file.                               *
	ParmParse pp("twocylinders");

	pp.query("direction1", direction1);
	pp.query("direction2", direction2);
	pp.query("radius1", radius1);
	pp.query("radius2", radius2);
	pp.getarr("center1", centervec1, 0, 3);
	pp.getarr("center2", centervec2, 0, 3);
	Array<Real, 3> center1 = {centervec1[0], centervec1[1], centervec1[2]};
	Array<Real, 3> center2 = {centervec2[0], centervec2[1], centervec2[2]};

    // Compute distance between cylinder centres
    Real offset = 0.0;
    for(int i = 0; i < 3; i++)
        offset += pow(center1[i] - center2[i], 2);
    offset = sqrt(offset); 

    // Print info about cylinders
	amrex::Print() << " CYLINDER 1" << std::endl;
	amrex::Print() << " Direction:       " << direction1 << std::endl;
	amrex::Print() << " Radius:    " << radius1 << std::endl;
    amrex::Print() << " Center:    " 
        << center1[0] << ", " << center1[1] << ", " << center1[2] << std::endl;

	amrex::Print() << " CYLINDER 2" << std::endl;
	amrex::Print() << " Direction:       " << direction2 << std::endl;
	amrex::Print() << " Radius:    " << radius2 << std::endl;
    amrex::Print() << " Center:    " 
        << center2[0] << ", " << center2[1] << ", " << center2[2] << std::endl;

	amrex::Print() << "\n Offset:          " << offset << std::endl;

	// Build the implicit function as a union of two cylinders
    EB2::CylinderIF cyl1(radius1, direction1, center1, false);
    EB2::CylinderIF cyl2(radius2, direction2, center2, false);
    auto twocylinders = EB2::makeUnion(cyl1, cyl2);

    // Generate GeometryShop
	auto gshop = EB2::makeShop(twocylinders);

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
