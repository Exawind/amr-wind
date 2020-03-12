#include <AMReX_EB2.H>
#include <AMReX_EB2_IF.H>

#include <AMReX_ParmParse.H>

#include <algorithm>
#include <incflo.H>

using namespace amrex;

/********************************************************************************
 *                                                                              *
 * Function to create two cylinders for the JCAP setup
 *                                                                              *
 ********************************************************************************/
void incflo::make_eb_cyl_tuscan()
{
    // Initialise parameters
    int direction1 = 0;
    int direction2 = 0;

    Real radius1 = 0.5;
    Real radius2 = 0.5;

    Real height1 = 0.3;
    Real height2 = 0.3;

    Vector<Real> centervec1(3);
    Vector<Real> centervec2(3);

    // Get information from inputs file.                               *
    ParmParse pp("jcap");

    pp.query("direction1", direction1);
    pp.query("direction2", direction2);

    pp.query("radius1", radius1);
    pp.query("radius2", radius2);

    pp.query("height1", height1);
    pp.query("height2", height2);

    pp.getarr("center1", centervec1, 0, 3);
    pp.getarr("center2", centervec2, 0, 3);

    Array<Real, 3> center1 = {centervec1[0], centervec1[1], centervec1[2]};
    Array<Real, 3> center2 = {centervec2[0], centervec2[1], centervec2[2]};

    // Print info about cylinders
    amrex::Print() << " CYLINDER 1" << std::endl;
    amrex::Print() << " Direction:       " << direction1 << std::endl;
    amrex::Print() << " Height:    " << height1 << std::endl;
    amrex::Print() << " Radius:    " << radius1 << std::endl;
    amrex::Print() << " Center:    " 
                   << center1[0] << ", " << center1[1] << ", " << center1[2] << std::endl;

    amrex::Print() << " CYLINDER 2" << std::endl;
    amrex::Print() << " Direction:       " << direction2 << std::endl;
    amrex::Print() << " Height:    " << height2 << std::endl;
    amrex::Print() << " Radius:    " << radius2 << std::endl;
    amrex::Print() << " Center:    " 
                   << center2[0] << ", " << center2[1] << ", " << center2[2] << std::endl;

    // Build the implicit function as a union of two cylinders
    EB2::CylinderIF cyl1(radius1, height1, direction1, center1, true);
    EB2::CylinderIF cyl2(radius2, height2, direction2, center2, true);
    auto twocylinders = EB2::makeIntersection(cyl1, cyl2);

    // Generate GeometryShop
    auto gshop = EB2::makeShop(twocylinders);

    // Build index space
    int max_level_here = 0;
    int max_coarsening_level = 100;
    EB2::Build(gshop, geom.back(), max_level_here, max_level_here + max_coarsening_level);
}
