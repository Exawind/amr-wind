#ifndef MOSD_H
#define MOSD_H

#include "amr-wind/boundary_conditions/wall_models/LogLaw.H"
#include "amr-wind/wind_energy/ShearStress.H"

namespace amr_wind {

struct MOSDSimple
{

    explicit MOSDSimple() {}

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
    get_shear(amrex::Real u, amrex::Real /* wspd */) const
    {
        0.;
    }
};

} // namespace amr_wind