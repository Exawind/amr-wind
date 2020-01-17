#include <AMReX_EB2.H>
#include <AMReX_EB2_IF.H>
#include <AMReX_ParmParse.H>

#include <algorithm>
#include <incflo.H>

using namespace amrex;

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
    EB2::Build(gshop, geom.back(), max_level_here, max_level_here + max_coarsening_level);
}
