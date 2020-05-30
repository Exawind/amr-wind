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
 *  Parameters are read from the `GeostrophicWind` and `CorolisForcing`
 *  namespace in the input file. The following parameters are available:
 *
 * - `rotational_time_period` Time period for planetary rotation (default: 86400
 *    seconds) in the CoriolisForcing namespace
 * 
 * - `geostrophic_wind` Geostrophic wind above capping inversion in the
 *    GeostrophicForcing namespace
 *
 */
GeostrophicForcing::GeostrophicForcing(const CFDSim& sim)
  : m_density(sim.repo().get_field("density"))
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

    m_g_forcing = {-coriolis_factor * m_target_vel[1],
                   coriolis_factor * m_target_vel[0],
                   0.0 };

}

GeostrophicForcing::~GeostrophicForcing() = default;

void GeostrophicForcing::operator()(
    const int lev,
    const amrex::MFIter& mfi,
    const amrex::Box& bx,
    const FieldState fstate,
    const amrex::Array4<amrex::Real>& src_term) const
{

    const auto& rho =
        m_density.state(fstate)(lev).const_array(mfi);
    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        src_term(i, j, k, 0) += rho(i,j,k) * m_g_forcing[0];
        src_term(i, j, k, 1) += rho(i,j,k) * m_g_forcing[1];
        // No forcing in z-direction
    });
}

} // namespace icns
} // namespace pde
} // namespace amr_wind
