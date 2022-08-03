#include "amr-wind/equation_systems/icns/source_terms/GeostrophicForcing.H"
#include "amr-wind/equation_systems/icns/source_terms/CoriolisForcing.H"
#include "amr-wind/CFDSim.H"
#include "amr-wind/utilities/trig_ops.H"
#include "amr-wind/core/vs/vstraits.H"

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
 *  - `rotational_time_period` Time period for planetary rotation (default:
 * 86400 seconds) in the CoriolisForcing namespace
 *
 *  - `geostrophic_wind` Geostrophic wind above capping inversion in the
 *    GeostrophicForcing namespace
 * 
 *  - 'three_ComponentForcing' (Default: false = 0 - two component forcing) 
 *
 */
GeostrophicForcing::GeostrophicForcing(const CFDSim& /*unused*/)
{
    amrex::Real coriolis_factor;
    {
        // Read the rotational time period (in seconds)
        amrex::ParmParse pp("CoriolisForcing");
        amrex::Real rot_time_period = 86400.0;
        pp.query("rotational_time_period", rot_time_period);
        coriolis_factor = 2.0 * utils::two_pi() / rot_time_period;
        amrex::Print() << "Geostrophic forcing: Coriolis factor = "
                       << coriolis_factor << std::endl;

        // Latitude is mandatory
        // Latitude is read in degrees
        pp.get("latitude", m_latitude);
        m_latitude = utils::radians(m_latitude);
        m_sinphi = std::sin(m_latitude);
        m_cosphi = std::cos(m_latitude);

        // 3-component forcing (Default: false)
        //bool m_S = false;
        if (!pp.query("three_ComponentForcing", m_S)){};
    }

    {
        // Read the geostrophic wind speed vector (in m/s)
        amrex::ParmParse pp("GeostrophicForcing");
        pp.getarr("geostrophic_wind", m_target_vel);
    }

    const auto sinphi = m_sinphi;
    const auto cosphi = m_cosphi;
    const auto corfac = m_coriolis_factor;
    const auto S = m_S;

    m_g_forcing = {
        -coriolis_factor * m_target_vel[1] * sinphi +coriolis_factor * m_target_vel[3] * cosphi * m_S, 
        +coriolis_factor * m_target_vel[0] * sinphi,
        -coriolis_factor * m_target_vel[0] * cosphi * m_S};
    amrex::Print() << "Geostrophic Forcing = " << m_g_forcing[0] <<", " << m_g_forcing[1] << ", " << m_g_forcing[2] << std::endl;
}

GeostrophicForcing::~GeostrophicForcing() = default;

void GeostrophicForcing::operator()(
    const int /*lev*/,
    const amrex::MFIter& /*mfi*/,
    const amrex::Box& bx,
    const FieldState /*fstate*/,
    const amrex::Array4<amrex::Real>& src_term) const
{
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> forcing{
        {m_g_forcing[0], m_g_forcing[1], m_g_forcing[2]}};
    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        src_term(i, j, k, 0) += forcing[0];
        src_term(i, j, k, 1) += forcing[1];
        src_term(i, j, k, 2) += forcing[2];
    });
}

} // namespace icns
} // namespace pde
} // namespace amr_wind
