#include "amr-wind/equation_systems/icns/source_terms/ABLForcing.H"
#include "amr-wind/CFDSim.H"
#include "amr-wind/wind_energy/ABL.H"
#include "amr-wind/physics/multiphase/MultiPhase.H"
#include "amr-wind/utilities/trig_ops.H"

#include "AMReX_ParmParse.H"
#include "AMReX_Gpu.H"

namespace amr_wind::pde::icns {

ABLForcing::ABLForcing(const CFDSim& sim) : m_time(sim.time()), m_mesh(sim.mesh())
{
    const auto& abl = sim.physics_manager().get<amr_wind::ABL>();
    abl.register_forcing_term(this);
    abl.abl_statistics().register_forcing_term(this);

    amrex::ParmParse pp_abl(identifier());
    // TODO: Allow forcing at multiple heights
    pp_abl.get("abl_forcing_height", m_forcing_height);
    amrex::ParmParse pp_incflo("incflo");

    pp_abl.query("velocity_timetable", m_vel_timetable);
    if (!m_vel_timetable.empty()) {
        std::ifstream ifh(m_vel_timetable, std::ios::in);
        if (!ifh.good()) {
            amrex::Abort("Cannot find input file: " + m_vel_timetable);
        }
        amrex::Real data_time;
        amrex::Real data_speed;
        amrex::Real data_deg;
        ifh.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        while (ifh >> data_time) {
            ifh >> data_speed >> data_deg;
            amrex::Real data_rad = utils::radians(data_deg);
            m_time_table.push_back(data_time);
            m_speed_table.push_back(data_speed);
            m_direction_table.push_back(data_rad);
        }
    } else {
        pp_incflo.getarr("velocity", m_target_vel);
    }

    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
        m_mean_vel[i] = m_target_vel[i];
    }

    // Set up relaxation toward 0 forcing near the air-water interface
    if (sim.repo().field_exists("vof")) {
        // If vof exists, get multiphase physics
        const auto& mphase = sim.physics_manager().get<amr_wind::MultiPhase>();
        // Retrieve interface position
        m_water_level = mphase.water_level();
        // Confirm that water level will be used
        m_use_phase_ramp = true;
        // Parse for thickness of ramping function
        pp_abl.get("abl_forcing_off_height", m_forcing_mphase0);
        pp_abl.get("abl_forcing_ramp_height", m_forcing_mphase1);
    }
}

ABLForcing::~ABLForcing() = default;

void ABLForcing::operator()(
    const int lev,
    const amrex::MFIter& /*mfi*/,
    const amrex::Box& bx,
    const FieldState /*fstate*/,
    const amrex::Array4<amrex::Real>& src_term) const
{

    const amrex::Real dudt = m_abl_forcing[0];
    const amrex::Real dvdt = m_abl_forcing[1];

    const bool ph_ramp = m_use_phase_ramp;
    const amrex::Real wlev = m_water_level;
    const amrex::Real wrht0 = m_forcing_mphase0;
    const amrex::Real wrht1 = m_forcing_mphase1;
    const auto& problo = m_mesh.Geom(lev).ProbLoArray();
    const auto& dx = m_mesh.Geom(lev).CellSizeArray();

    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        amrex::Real fac = 1.0;
        if (ph_ramp) {
            const amrex::Real z = problo[2] + (k + 0.5) * dx[2];
            if (z - wlev < wrht0 + wrht1) {
                if (z - wlev < wrht0) {
                    // Apply no forcing within first interval
                    fac = 0.0;
                } else {
                    // Ramp from 0 to 1 over second interval
                    fac =
                        0.5 - 0.5 * std::cos(M_PI * (z - wlev - wrht0) / wrht1);
                }
            }
        }
        src_term(i, j, k, 0) += fac * dudt;
        src_term(i, j, k, 1) += fac * dvdt;

        // No forcing in z-direction
    });
}

} // namespace amr_wind::pde::icns
