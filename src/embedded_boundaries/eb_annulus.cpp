#include <AMReX_EB2.H>
#include <AMReX_EB2_IF_Cylinder.H>
#include <AMReX_EB2_IF_Plane.H>
#include <AMReX_EB2_IF_Union.H>

#include <AMReX_ParmParse.H>

#include <AMReX_EB_levelset.H>
#include <algorithm>
#include <embedded_boundaries_F.H>
#include <incflo.H>

/********************************************************************************
 *                                                                              *
 * Function to create a annular cylinder EB.                                     *
 *                                                                              *
 ********************************************************************************/
void incflo::make_eb_annulus()
{
	ParmParse pp("annulus");

	int max_level_here = 0;

	/****************************************************************************
     * Get cylinder information from inputs file.                               *
     ****************************************************************************/
	int direction = 0;
	Real height = -1.0;
	Real outer_radius = 0.0002;
	Real inner_radius = 0.0001;
	Vector<Real> outer_centervec(3);
	Vector<Real> inner_centervec(3);

	pp.query("direction", direction);
	pp.query("height", height);
	pp.query("outer_radius", outer_radius);
	pp.query("inner_radius", inner_radius);
	pp.getarr("outer_center", outer_centervec, 0, 3);
	pp.getarr("inner_center", inner_centervec, 0, 3);
	Array<Real, 3> outer_center = {outer_centervec[0], outer_centervec[1], outer_centervec[2]};
	Array<Real, 3> inner_center = {inner_centervec[0], inner_centervec[1], inner_centervec[2]};

    // make_eb_annulus: outer and inner cylinders must have the same center coordinate in given direction
    AMREX_ASSERT(outer_center[direction] == inner_center[direction]);

    // Compute distance between cylinder centres
    Real offset = 0.0;
    for(int i = 0; i < 3; i++)
        offset += pow(outer_center[i] - inner_center[i], 2);
    offset = sqrt(offset); 

    // Check that the inner cylinder is fully contained in the outer one
    Real smallest_gap_width = outer_radius - inner_radius - offset;
    AMREX_ASSERT(smallest_gap_width >= 0.0);

    // Compute standoff - measure of eccentricity
    Real standoff = 100 * smallest_gap_width / (outer_radius - inner_radius);
    AMREX_ASSERT((standoff >= 0) && (standoff <= 100));
    

	/****************************************************************************
     *                                                                          *
     * Build standard EB Factories                                              *
     *                                                                          *
     ****************************************************************************/

	// set up ebfactory

	EBSupport m_eb_support_level = EBSupport::full;

	amrex::Print() << " " << std::endl;
	amrex::Print() << " Direction:       " << direction << std::endl;
	amrex::Print() << " Height:          " << ((height < 0.0) ? "inf" : std::to_string(height)) << std::endl;
	amrex::Print() << " Outer radius:    " << outer_radius << std::endl;
	amrex::Print() << " Inner radius:    " << inner_radius << std::endl;
    amrex::Print() << " Outer center:    " 
        << outer_center[0] << ", " << outer_center[1] << ", " << outer_center[2] << std::endl;
    amrex::Print() << " Inner center:    " 
        << inner_center[0] << ", " << inner_center[1] << ", " << inner_center[2] << std::endl;
	amrex::Print() << " Offset:          " << offset << std::endl;
	amrex::Print() << " Smallest gap:    " << smallest_gap_width << std::endl;
	amrex::Print() << " Standoff:        " << standoff << std::endl;

	// Create the annulus
	amrex::Print() << "Building the annular cylinder geometry ..." << std::endl;

    std::unique_ptr<EB2::CylinderIF> outer_cyl;
    std::unique_ptr<EB2::CylinderIF> inner_cyl;
    if (height < 0.0) 
    {
        outer_cyl = std::unique_ptr<EB2::CylinderIF>(
                new EB2::CylinderIF(outer_radius, direction, outer_center, true));
        inner_cyl = std::unique_ptr<EB2::CylinderIF>(
                new EB2::CylinderIF(inner_radius, direction, inner_center, false));
    } 
    else 
    {
        outer_cyl = std::unique_ptr<EB2::CylinderIF>(
                new EB2::CylinderIF(outer_radius, height, direction, outer_center, true));
        inner_cyl = std::unique_ptr<EB2::CylinderIF>(
                new EB2::CylinderIF(inner_radius, height, direction, inner_center, false));
    }

    auto annulus = EB2::makeUnion(*outer_cyl, *inner_cyl);

	auto gshop = EB2::makeShop(annulus);
	int max_coarsening_level = 100;
	EB2::Build(gshop, geom.back(), max_level_here, max_level_here + max_coarsening_level);

    for(int lev = 0; lev < nlev; lev++)
    {
        const EB2::IndexSpace& ebis = EB2::IndexSpace::top();
        const EB2::Level& ebis_lev = ebis.getLevel(geom[lev]);

        amrex::Print() << "Done building the annular cylinder geometry" << std::endl;

        /****************************************************************************
        *                                                                           *
        * THIS FILLS FLUID EBFACTORY                                                *
        *                                                                           *
        *****************************************************************************/

        amrex::Print() << "Now  making the fluid ebfactory ..." << std::endl;

        eb_level_fluid = &ebis_lev;

        ebfactory[lev].reset(new EBFArrayBoxFactory(
            *eb_level_fluid,
            geom[lev],
            grids[lev],
            dmap[lev],
            {m_eb_basic_grow_cells, m_eb_volume_grow_cells, m_eb_full_grow_cells},
            m_eb_support_level));

        amrex::Print() << "Done making the fluid ebfactory ..." << std::endl;
    }
}
