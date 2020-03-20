#include "ABLWallFunction.H"

#include <cmath>

#include "AMReX_ParmParse.H"
#include "AMReX_Print.H"

namespace amr_wind {

ABLWallFunction::ABLWallFunction()
{
    amrex::ParmParse pp("abl");

    pp.query("kappa", m_kappa);
    pp.query("surface_roughness_z0", m_z0);
}

amrex::Real ABLWallFunction::utau(const amrex::Real uvel, const amrex::Real zh) const
{
    return m_kappa * uvel / std::log(zh / m_z0);
}

}
