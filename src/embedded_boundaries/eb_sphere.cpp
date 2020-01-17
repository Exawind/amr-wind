#include <AMReX_EB2.H>
#include <AMReX_EB2_IF.H>
#include <AMReX_ParmParse.H>

#include <algorithm>
#include <incflo.H>

using namespace amrex;

/********************************************************************************
 *                                                                              *
 * Function to create a simple sphere EB.                                     *
 *                                                                              *
 ********************************************************************************/
void incflo::make_eb_sphere()
{
    // Initialise sphere parameters
    bool inside = true;
    Real radius = 0.0002;
    Vector<Real> centervec(3);

    // Get sphere information from inputs file.                               *
    ParmParse pp("sphere");

    pp.query("internal_flow", inside);
    pp.query("radius", radius);
    pp.getarr("center", centervec, 0, 3);
    Array<Real, 3> center = {centervec[0], centervec[1], centervec[2]};

    // Print info about sphere
    amrex::Print() << " " << std::endl;
    amrex::Print() << " Internal Flow: " << inside << std::endl;
    amrex::Print() << " Radius:    " << radius << std::endl;
    amrex::Print() << " Center:    " << center[0] << ", " << center[1] << ", " << center[2]
                   << std::endl;

    // Build the sphere implicit function 
    EB2::SphereIF my_sphere(radius, center, inside);

    // Generate GeometryShop
    auto gshop = EB2::makeShop(my_sphere);

    // Build index space
    int max_level_here = 0;
    int max_coarsening_level = 100;
    EB2::Build(gshop, geom.back(), max_level_here, max_level_here + max_coarsening_level);
}
