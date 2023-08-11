#include "amr-wind/equation_systems/icns/source_terms/GeostrophicForcing.H"
#include "amr-wind/CFDSim.H"
#include "amr-wind/utilities/trig_ops.H"
#include "amr-wind/core/vs/vstraits.H"

#include "AMReX_ParmParse.H"
#include "AMReX_Gpu.H"

namespace amr_wind::pde::icns {

/** Geostrophic forcing term for ABL
 *
 *  Parameters are read from the `GeostrophicWind` and `CorolisForcing`
 *  namespace in the input file. The following parameters are available:
 *
 *  - `rotational_time_period` Time period for planetary rotation (default:
 * 86400 seconds) in the CoriolisForcing namespace
 *
 *  - `geostrophic_wind` Geostrophic wind above capping inversion in the
 *    GeostrophicForcing namespace
 *
 */
GeostrophicForcing::GeostrophicForcing(const CFDSim& /*unused*/)
{
    amrex::Real coriolis_factor = 0.0;

    // Read the rotational time period (in seconds)
    amrex::ParmParse ppc("CoriolisForcing");
    // Read the rotational time period (in seconds) -- This is 23hrs and 56
    // minutes and 4.091 seconds
    amrex::Real rot_time_period = 86164.091;
    ppc.query("rotational_time_period", rot_time_period);

    amrex::Real omega = 2.0 * utils::pi() / rot_time_period;
    amrex::Real latitude = 90;
    ppc.get("latitude", latitude);
    latitude = utils::radians(latitude);
    amrex::Real sinphi = std::sin(latitude);

    coriolis_factor = 2.0 * omega * sinphi;
    ppc.query("is_horizontal", m_is_horizontal);
    amrex::Print() << "Geostrophic forcing: Coriolis factor = "
                   << coriolis_factor << std::endl;

    // Read the geostrophic wind speed vector (in m/s)
    amrex::ParmParse ppg("GeostrophicForcing");
    ppg.getarr("geostrophic_wind", m_target_vel);

    m_g_forcing = {
        -coriolis_factor * m_target_vel[1], coriolis_factor * m_target_vel[0],
        0.0};
}

GeostrophicForcing::~GeostrophicForcing() = default;

void GeostrophicForcing::operator()(
    const int /*lev*/,
    const amrex::MFIter& /*mfi*/,
    const amrex::Box& bx,
    const FieldState /*fstate*/,
    const amrex::Array4<amrex::Real>& src_term) const
{
    amrex::Real fac = (m_is_horizontal) ? 0. : 1.;
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> forcing{
        {m_g_forcing[0], m_g_forcing[1], m_g_forcing[2]}};
    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        src_term(i, j, k, 0) += forcing[0];
        src_term(i, j, k, 1) += forcing[1];
        src_term(i, j, k, 1) += fac * forcing[2];
    });
}

} // namespace amr_wind::pde::icns
