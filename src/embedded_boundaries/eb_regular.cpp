#include <AMReX_EB2.H>
#include <AMReX_EB2_IF_Plane.H>
#include <AMReX_EB2_IF_Union.H>

#include <eb_if.H> 
#include <incflo.H>

/********************************************************************************
 *                                                                              *
 * Placeholder: create a simulation box _without_ EB walls.                     *
 *                                                                              *
 ********************************************************************************/
void incflo::make_eb_regular()
{
    // Check whether there are walls present
    bool has_walls = false;
    // std::unique_ptr<UnionListIF<EB2::PlaneIF>> impfunc_walls = get_real_walls(has_walls);

    // Settings for building index space
    int max_level_here = 0;
    int max_coarsening_level = 100;

    // Generate GeometryShop and build index space
#if 0
    if(has_walls)
    {
        auto gshop = EB2::makeShop(*impfunc_walls);
        EB2::Build(gshop, geom.back(), max_level_here, max_level_here + max_coarsening_level);
    }
    else
#endif
    {
        EB2::AllRegularIF my_regular;
        auto gshop = EB2::makeShop(my_regular);
        EB2::Build(gshop, geom.back(), max_level_here, max_level_here + max_coarsening_level);
    }
}
