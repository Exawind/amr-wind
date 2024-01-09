#include "amr-wind/equation_systems/icns/source_terms/GeostrophicForcing.H"
#include "amr-wind/CFDSim.H"
#include "amr-wind/utilities/trig_ops.H"
#include "amr-wind/core/vs/vstraits.H"
#include "amr-wind/physics/multiphase/MultiPhase.H"
#include "amr-wind/equation_systems/vof/volume_fractions.H"

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
GeostrophicForcing::GeostrophicForcing(const CFDSim& sim) : m_mesh(sim.mesh())
{
    amrex::Real coriolis_factor = 0.0;

    // Read the rotational time period (in seconds)
    amrex::ParmParse ppc("CoriolisForcing");
    // Read the rotational time period (in seconds) -- This is 23hrs and 56
    // minutes and 4.091 seconds
    amrex::Real rot_time_period = 86164.091;
    ppc.query("rotational_time_period", rot_time_period);

    amrex::Real omega = utils::two_pi() / rot_time_period;
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

    // Set up relaxation toward 0 forcing near the air-water interface
    if (sim.repo().field_exists("vof")) {
        // If vof exists, get multiphase physics
        const auto& mphase = sim.physics_manager().get<amr_wind::MultiPhase>();
        // Retrieve interface position
        m_water_level = mphase.water_level();
        // Confirm that water level will be used
        m_use_phase_ramp = true;
        // Parse for thickness of ramping function
        ppg.get("wind_forcing_off_height", m_forcing_mphase0);
        ppg.get("wind_forcing_ramp_height", m_forcing_mphase1);
        // Store reference to vof field
        m_vof = &sim.repo().get_field("vof");
        // Parse for number of cells in band
        ppg.query("wind_forcing_band", m_n_band);
    }
}

GeostrophicForcing::~GeostrophicForcing() = default;

void GeostrophicForcing::operator()(
    const int lev,
    const amrex::MFIter& mfi,
    const amrex::Box& bx,
    const FieldState /*fstate*/,
    const amrex::Array4<amrex::Real>& src_term) const
{
    amrex::Real hfac = (m_is_horizontal) ? 0. : 1.;

    const bool ph_ramp = m_use_phase_ramp;
    const int n_band = m_n_band;
    const amrex::Real wlev = m_water_level;
    const amrex::Real wrht0 = m_forcing_mphase0;
    const amrex::Real wrht1 = m_forcing_mphase1;
    const auto& problo = m_mesh.Geom(lev).ProbLoArray();
    const auto& dx = m_mesh.Geom(lev).CellSizeArray();
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> forcing{
        {m_g_forcing[0], m_g_forcing[1], m_g_forcing[2]}};

    const auto& vof = (*m_vof)(lev).const_array(mfi);
    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        amrex::Real wfac = 1.0;
        if (ph_ramp) {
            const amrex::Real z = problo[2] + (k + 0.5) * dx[2];
            if (z - wlev < wrht0 + wrht1) {
                if (z - wlev < wrht0) {
                    // Apply no forcing within first interval
                    wfac = 0.0;
                } else {
                    // Ramp from 0 to 1 over second interval
                    wfac =
                        0.5 - 0.5 * std::cos(M_PI * (z - wlev - wrht0) / wrht1);
                }
            }
            // Check for presence of liquid (like a droplet)
            // - interface_band checks for closeness to interface
            // - need to also check for within liquid
            if (multiphase::interface_band(i, j, k, vof, n_band) ||
                vof(i, j, k) > 1.0 - 1e-12) {
                // Turn off forcing
                wfac = 0.0;
            }
        }
        src_term(i, j, k, 0) += wfac * forcing[0];
        src_term(i, j, k, 1) += wfac * forcing[1];
        src_term(i, j, k, 1) += wfac * hfac * forcing[2];
    });
}

} // namespace amr_wind::pde::icns
