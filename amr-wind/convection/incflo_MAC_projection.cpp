#include <AMReX_REAL.H>
#include <AMReX_BLFort.H>
#include <AMReX_SPACE.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_MacProjector.H>
#include <AMReX_Vector.H>

#include "amr-wind/convection/mac_projection.H"
#include "amr-wind/core/MLMGOptions.H"
#include "amr-wind/utilities/console_io.H"

using namespace amrex;

Array<amrex::LinOpBCType,AMREX_SPACEDIM>
mac::get_projection_bc (Orientation::Side side, GpuArray<BC, AMREX_SPACEDIM*2> bctype, Vector<Geometry> geom) noexcept
{

    Array<LinOpBCType,AMREX_SPACEDIM> r;
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        if (geom[0].isPeriodic(dir)) {
            r[dir] = LinOpBCType::Periodic;
        } else {
            auto bc = bctype[Orientation(dir,side)];
            switch (bc)
            {
            case BC::pressure_inflow:
            case BC::pressure_outflow:
            {
                r[dir] = LinOpBCType::Dirichlet;
                break;
            }
            default:
                r[dir] = LinOpBCType::Neumann;
                break;
            };
        }
    }
    return r;
}

