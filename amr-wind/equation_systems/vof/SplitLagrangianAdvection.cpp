#include "amr-wind/equation_systems/vof/SplitLagrangianAdvection.H"
#include <AMReX_Geometry.H>

using namespace amrex;

namespace amr_wind {

void multiphase::split_lagrangian_advection(
    int lev,
    amrex::Box const& bx,
    amrex::Array4<amrex::Real const> const& q,
    amrex::Array4<amrex::Real const> const& umac,
    amrex::Array4<amrex::Real const> const& vmac,
    amrex::Array4<amrex::Real const> const& wmac,
    amrex::BCRec const* pbc,
    amrex::Vector<amrex::Geometry> geom,
    amrex::Real dt)
{
    BL_PROFILE("amr-wind::godunov::compute_advection");
    Box const& xbx = amrex::surroundingNodes(bx, 0);
    Box const& ybx = amrex::surroundingNodes(bx, 1);
    Box const& zbx = amrex::surroundingNodes(bx, 2);
    Box const& bxg1 = amrex::grow(bx, 1);
    Box xebox = Box(xbx).grow(1, 1).grow(2, 1);
    Box yebox = Box(ybx).grow(0, 1).grow(2, 1);
    Box zebox = Box(zbx).grow(0, 1).grow(1, 1);

    const Real dx = geom[lev].CellSize(0);
    const Real dy = geom[lev].CellSize(1);
    const Real dz = geom[lev].CellSize(2);
    Real l_dt = dt;
    Real dtdx = l_dt / dx;
    Real dtdy = l_dt / dy;
    Real dtdz = l_dt / dz;
}

} // namespace amr_wind