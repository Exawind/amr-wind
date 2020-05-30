#include "GeostrophicForcing.H"
#include "PlaneAveraging.H"
#include "CFDSim.H"
#include "trig_ops.H"

#include "AMReX_ParmParse.H"
#include "AMReX_Gpu.H"

namespace amr_wind {
namespace pde {
namespace icns {

/** Geostrophic forcing term for ABL 
 *
 *  Parameters are read from the `GeostrophicWind`, `CorolisForcing`, and
 *  `incflo` namespace in the input file. The following parameters are
 *  available:
 *
 * - `rotational_time_period` Time period for planetary rotation (default: 86400
 *    seconds) in the CoriolisForcing namespace
 * 
 * - `geostrophic_wind` Geostrophic wind above capping inversion in the
 *    GeostrophicForcing namespace
 *
 * - `density` Density in the incflo namespace
 *
 */
GeostrophicForcing::GeostrophicForcing(const CFDSim&)
{
    amrex::Real coriolis_factor;
    {
        // Read the rotational time period (in seconds)
        amrex::ParmParse pp("CoriolisForcing");
        amrex::Real rot_time_period = 86400.0;
        pp.query("rotational_time_period", rot_time_period);
        coriolis_factor = 2.0 * utils::two_pi() / rot_time_period;
    }
    
    {
        // Read the geostrophic wind speed vector (in m/s)
        amrex::ParmParse pp("GeostrophicForcing");
        pp.getarr("geostrophic_wind", m_target_vel);
    }

    amrex::Real rho=1.0;
    {
        // Read the density
        amrex::ParmParse pp("incflo");
        pp.query("density", rho);
    }
    
    m_g_forcing = {-rho* coriolis_factor * m_target_vel[1],
                   rho * coriolis_factor * m_target_vel[0],
                   0.0 };

}

GeostrophicForcing::~GeostrophicForcing() = default;

void GeostrophicForcing::operator()(
    const int,
    const amrex::MFIter&,
    const amrex::Box& bx,
    const FieldState,
    const amrex::Array4<amrex::Real>& src_term) const
{

    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        src_term(i, j, k, 0) += m_g_forcing[0];
        src_term(i, j, k, 1) += m_g_forcing[1];
        // No forcing in z-direction
    });
}

} // namespace icns
} // namespace pde
} // namespace amr_wind
