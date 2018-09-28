//#include <AMReX_EB2.H>
//#include <AMReX_EB2_IF_Cylinder.H>
//#include <AMReX_EB2_IF_Plane.H>
//#include <AMReX_EB2_IF_Union.H>

//#include <AMReX_VisMF.H>  // amrex::VisMF::Write(MultiFab)
//#include <sstream>

#include <AMReX_ParmParse.H>

#include <AMReX_EB_levelset.H>
#include <algorithm>
#include <eb_F.H>
#include <incflo_level.H>

void incflo_level::make_eb_geometry(int lev)
{

	/******************************************************************************
   * incflo.geometry=<string> specifies the EB geometry. <string> can be on of    *
   * box, cylinder, hopper, clr, clr_riser, general (or blank)                  *
   ******************************************************************************/

	ParmParse pp("incflo");

	std::string geom_type;
	pp.query("geometry", geom_type);

	/******************************************************************************
   * Legacy inputs:                                                             *
   *   -- incflo.hourglass = true <=> incflo.geometry=box                           *
   *   -- incflo.clr       = true <=> incflo.geometry=clr                           *
   *   -- incflo.clr_riser = true <=> incflo.geometry=clr_riser                     *
   *   -- incflo.use_walls = true <=> incflo.geometry=general                       *
   *   -- incflo.use_poy2  = true <=> incflo.geometry=general                       *
   ******************************************************************************/

	bool hourglass = false;
	bool clr = false;
	bool clr_riser = false;
	bool eb_general = false;

	pp.query("hourglass", hourglass);
	pp.query("clr", clr);
	pp.query("clr_riser", clr_riser);

	bool eb_poly2 = false;
	bool eb_walls = false;

	pp.query("use_poly2", eb_poly2);
	pp.query("use_walls", eb_walls);
	eb_general = eb_poly2 || eb_walls;

	// Avoid multiple (ambiguous) inputs
	if(hourglass || clr || clr_riser || eb_general)
	{
		if(!geom_type.empty())
		{
			amrex::Abort("The input file cannot specify both:\n"
						 "incflo.<geom_type>=true and incflo.geometry=<geom_type>\n"
						 "at the same time.");
		}
	}

	if(hourglass)
		geom_type = "hourglass";
	if(clr)
		geom_type = "clr";
	if(clr_riser)
		geom_type = "clr_riser";
	if(eb_general)
		geom_type = "general";

	/******************************************************************************
   *                                                                            *
   *  CONSTRUCT EB                                                              *
   *                                                                            *
   ******************************************************************************/

	if(geom_type == "box")
	{
		amrex::Print() << "\n Building box geometry." << std::endl;
		make_eb_box(lev);
	}
	else if(geom_type == "annulus")
	{
		amrex::Print() << "\n Building annulus geometry." << std::endl;
		make_eb_annulus(lev);
	}
	else if(geom_type == "cylinder")
	{
		amrex::Print() << "\n Building cylinder geometry." << std::endl;
		make_eb_cylinder(lev);
	}
	else if(geom_type == "hopper")
	{
		amrex::Print() << "\n Building hopper geometry." << std::endl;
		make_eb_hopper(lev);
	}
	else if(geom_type == "general")
	{
		amrex::Print() << "\n Building general geometry (poly2 with extra walls)." << std::endl;
		make_eb_general(lev);
	}
	else
	{
		amrex::Print() << "\n No EB geometry declared in inputs => "
					   << " Will read walls from incflo.dat only." << std::endl;
		make_eb_regular(lev);
	}
}
