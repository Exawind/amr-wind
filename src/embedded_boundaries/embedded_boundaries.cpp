#include <AMReX_ParmParse.H>

#include <algorithm>
#include <incflo.H>

using namespace amrex;

void incflo::MakeEBGeometry()
{
   /******************************************************************************
   * incflo.geometry=<string> specifies the EB geometry. <string> can be one of    *
   * box, cylinder, annulus, sphere, spherecube, twocylinders
   ******************************************************************************/

    ParmParse pp("incflo");

    std::string geom_type;
    pp.query("geometry", geom_type);

   /******************************************************************************
   *                                                                            *
   *  CONSTRUCT EB                                                              *
   *                                                                            *
   ******************************************************************************/

    if(geom_type == "box")
    {
	amrex::Print() << "\n Building box geometry." << std::endl;
        make_eb_box();
    }
    else if(geom_type == "cylinder")
    {
	amrex::Print() << "\n Building cylinder geometry." << std::endl;
        make_eb_cylinder();
    }
    else if(geom_type == "twocylinders")
    {
	amrex::Print() << "\n Building twocylinders geometry." << std::endl;
        make_eb_twocylinders();
    }
    else if(geom_type == "annulus")
    {
	amrex::Print() << "\n Building annulus geometry." << std::endl;
        make_eb_annulus();
    }
    else if(geom_type == "sphere")
    {
	amrex::Print() << "\n Building sphere geometry." << std::endl;
        make_eb_sphere();
    }
    else if(geom_type == "spherecube")
    {
	amrex::Print() << "\n Building spherecube geometry." << std::endl;
        make_eb_spherecube();
    }
    else if(geom_type == "tuscan")
    {
	amrex::Print() << "\n Building tuscan geometry." << std::endl;
        make_eb_tuscan();
    }
    else if(geom_type == "jcap")
    {
	amrex::Print() << "\n Building JCAP geometry." << std::endl;
        make_eb_cyl_tuscan();
    }
    else
    {
	amrex::Print() << "\n No EB geometry declared in inputs => "
	               << " Will build all regular geometry." << std::endl;
        make_eb_regular();
    }
    amrex::Print() << "Done making the geometry ebfactory.\n" << std::endl;
}
