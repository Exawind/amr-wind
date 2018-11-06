#include <AMReX_EB2.H>
#include <AMReX_EB2_IF_Cylinder.H>
#include <AMReX_EB2_IF_Plane.H>
#include <AMReX_EB2_IF_Union.H>

#include <algorithm>
#include <eb_if.H>
#include <embedded_boundaries_F.H>
#include <incflo.H>

std::unique_ptr<UnionListIF<EB2::PlaneIF>> incflo::get_walls(bool& has_walls)
{
	// Extracts all walls from the incflo.dat

	has_walls = false; // will be set to true if there are any walls

	// Walls can be defined per phase => Itterarte over all phases and check
	// each for walls in the incflo.dat
	Vector<EB2::PlaneIF> planes;
	for(int i = 1; i <= 500; i++)
	{
		int exists;
		RealVect normal, center;
		incflo_get_walls(&i, &exists, &normal, &center);
		if(exists)
		{
			has_walls = true;
			amrex::Print() << "Normal " << normal << std::endl;
			amrex::Print() << "Center " << center << std::endl;

			RealArray plane_center = {AMREX_D_DECL(center[0], center[1], center[2])};
			RealArray plane_normal = {AMREX_D_DECL(normal[0], normal[1], normal[2])};

			planes.emplace_back(plane_center, plane_normal, false);
		}
	}

	std::unique_ptr<UnionListIF<EB2::PlaneIF>> ret =
		std::unique_ptr<UnionListIF<EB2::PlaneIF>>(new UnionListIF<EB2::PlaneIF>(planes));
	return ret;
}

std::unique_ptr<UnionListIF<EB2::PlaneIF>> incflo::get_real_walls(bool& has_real_walls)
{
	// Extracts all walls from the incflo.dat

	has_real_walls = false; // will be set to true if there are any walls

	// Walls can be defined per phase => Itterarte over all phases and check
	// each for walls in the incflo.dat
	Vector<EB2::PlaneIF> planes;
	for(int i = 1; i <= 500; i++)
	{
		int exists;
		RealVect normal, center;
		incflo_get_real_walls(&i, &exists, &normal, &center);
		if(exists)
		{
			has_real_walls = true;
			amrex::Print() << "Normal " << normal << std::endl;
			amrex::Print() << "Center " << center << std::endl;

			RealArray plane_center = {AMREX_D_DECL(center[0], center[1], center[2])};
			RealArray plane_normal = {AMREX_D_DECL(normal[0], normal[1], normal[2])};

			planes.emplace_back(plane_center, plane_normal, false);
		}
	}

	std::unique_ptr<UnionListIF<EB2::PlaneIF>> ret =
		std::unique_ptr<UnionListIF<EB2::PlaneIF>>(new UnionListIF<EB2::PlaneIF>(planes));
	return ret;
}
