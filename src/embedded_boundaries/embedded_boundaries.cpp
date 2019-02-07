#include <AMReX_ParmParse.H>

#include <algorithm>
#include <embedded_boundaries_F.H>
#include <incflo.H>

void incflo::MakeEBGeometry()
{
    MakeBCArrays();

	/******************************************************************************
   * incflo.geometry=<string> specifies the EB geometry. <string> can be one of    *
   * box, cylinder, sphere, general (or blank)                                     *
   ******************************************************************************/

	ParmParse pp("incflo");

	std::string geom_type;
	pp.query("geometry", geom_type);

	/******************************************************************************
   * Legacy inputs:                                                             *
   *   -- incflo.use_walls = true <=> incflo.geometry=general                       *
   *   -- incflo.use_poy2  = true <=> incflo.geometry=general                       *
   ******************************************************************************/

	bool eb_general = false;

	bool eb_poly2 = false;
	bool eb_walls = false;

	pp.query("use_poly2", eb_poly2);
	pp.query("use_walls", eb_walls);
	eb_general = eb_poly2 || eb_walls;

	// Avoid multiple (ambiguous) inputs
	if(eb_general)
	{
		if(!geom_type.empty())
		{
			amrex::Abort("The input file cannot specify both:\n"
						 "incflo.<geom_type>=true and incflo.geometry=<geom_type>\n"
						 "at the same time.");
		}
	}

    if (eb_general) geom_type = "general";

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
	else if(geom_type == "annulus")
	{
		amrex::Print() << "\n Building annulus geometry." << std::endl;
        make_eb_annulus();
	}
	else if(geom_type == "cylinder")
	{
		amrex::Print() << "\n Building cylinder geometry." << std::endl;
        make_eb_cylinder();
	}
	else if(geom_type == "sphere")
	{
		amrex::Print() << "\n Building sphere geometry." << std::endl;
        make_eb_sphere();
	}
	else if(geom_type == "general")
	{
		amrex::Print() << "\n Building general geometry (poly2 with extra walls)." << std::endl;
        make_eb_general();
	}
	else
	{
		amrex::Print() << "\n No EB geometry declared in inputs => "
					   << " Will read walls from incflo.dat only." << std::endl;
        make_eb_regular();
	}
}

// This function checks if ebfactory is allocated with
// the proper dm and ba
bool incflo::UpdateEBFactory(int a_lev)
{
    // This assert is to verify that some kind of EB geometry
    // has already been defined
    AMREX_ASSERT(not EB2::IndexSpace::empty());

    const DistributionMapping& dm = DistributionMap(a_lev);
    const BoxArray&            ba = boxArray(a_lev);
    const EB2::IndexSpace&   ebis = EB2::IndexSpace::top();
    const EB2::Level&  ebis_level = ebis.getLevel(geom[a_lev]);

    bool is_updated = false;

    EBSupport m_eb_support_level = EBSupport::full;
    if ( ebfactory[a_lev].get() == nullptr )
    {
        ebfactory[a_lev].reset(new EBFArrayBoxFactory(ebis_level, geom[a_lev], ba, dm,
                                                      {m_eb_basic_grow_cells,
                                                       m_eb_volume_grow_cells,
                                                       m_eb_full_grow_cells},
                                                       m_eb_support_level));
        is_updated = true;
    }
    else
    {
        const DistributionMapping&  eb_dm = ebfactory[a_lev]->DistributionMap();
        const BoxArray&             eb_ba = ebfactory[a_lev]->boxArray();

        if ( (dm != eb_dm) || (ba != eb_ba) )
        {

            ebfactory[a_lev].reset(new EBFArrayBoxFactory(ebis_level, geom[a_lev], ba, dm,
                                                          {m_eb_basic_grow_cells,
                                                           m_eb_volume_grow_cells,
                                                           m_eb_full_grow_cells}, 
                                                           m_eb_support_level));
            is_updated = true;
        }
    }

    return is_updated;
}
