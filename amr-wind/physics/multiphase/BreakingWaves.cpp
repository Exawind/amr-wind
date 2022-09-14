#include "amr-wind/physics/multiphase/MultiPhase.H"
#include "amr-wind/physics/multiphase/BreakingWaves.H"
#include "amr-wind/utilities/trig_ops.H"
#include "amr-wind/CFDSim.H"
#include "AMReX_ParmParse.H"

namespace amr_wind {

BreakingWaves::BreakingWaves(CFDSim& sim)
    : m_sim(sim)
    , m_velocity(sim.repo().get_field("velocity"))
    , m_levelset(sim.repo().get_field("levelset"))
    , m_density(sim.repo().get_field("density"))
{

    {
        std::string turbulence_model;
        amrex::ParmParse pp("turbulence");
        pp.query("model", turbulence_model);
        if (turbulence_model == "Laminar") {
            m_laminar = true;
        }
    }

    {
        amrex::ParmParse pp("BreakingWaves");
        pp.query("amplitude", m_amplitude);
        pp.query("wavelength", m_wavelength);
        pp.query("water_surface_level", m_waterlevel);

        if (m_laminar) {
            pp.query("air_vel_shear_rate", m_air_vel_shear_rate);
        } else {
            pp.query("air_vel_shear_rate", m_air_vel_shear_rate);
            pp.query("tke0", m_tke0);
            pp.query("sdr0", m_sdr0);
        }
    }
}

BreakingWaves::~BreakingWaves() = default;

/** Initialize the velocity and vof fields at the beginning of the
 *  simulation.
 */
void BreakingWaves::initialize_fields(int level, const amrex::Geometry& geom)
{

    auto& velocity = m_velocity(level);
    auto& levelset = m_levelset(level);
    auto& density = m_density(level);

    const auto& mphase = m_sim.physics_manager().get<MultiPhase>();
    const amrex::Real rho1 = mphase.rho1();
    const amrex::Real rho2 = mphase.rho2();
    velocity.setVal(0.0, 0, AMREX_SPACEDIM);

    const auto& dx = geom.CellSizeArray();
    const auto& problo = geom.ProbLoArray();
    const amrex::Real alpha = m_amplitude;
    const amrex::Real lambda = m_wavelength;
    const amrex::Real water_level = m_waterlevel;
    const amrex::Real vel_shear_rate = m_air_vel_shear_rate;

    if (!m_laminar) {
        auto& tke = m_sim.repo().get_field("tke")(level);
        auto& sdr = m_sim.repo().get_field("sdr")(level);
        tke.setVal(m_tke0);
        sdr.setVal(m_sdr0);
    }

    for (amrex::MFIter mfi(levelset); mfi.isValid(); ++mfi) {
        const auto& vbx = mfi.growntilebox();
        auto vel = velocity.array(mfi);
        auto phi = levelset.array(mfi);
        auto rho = density.array(mfi);

        amrex::ParallelFor(
            vbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                const amrex::Real x = problo[0] + (i + 0.5) * dx[0];
                const amrex::Real z = problo[2] + (k + 0.5) * dx[2];
                const amrex::Real kappa = 2.0 * utils::pi() / lambda;
                const amrex::Real epsilon = alpha * kappa;
                // Compute free surface amplitude
                const amrex::Real eta =
                    water_level +
                    alpha * ((1.0 - 1.0 / 16.0 * epsilon * epsilon) *
                                 std::cos(kappa * x) +
                             0.5 * epsilon * std::cos(2.0 * kappa * x) +
                             3.0 / 8.0 * epsilon * epsilon *
                                 std::cos(3.0 * kappa * x));
                phi(i, j, k) = eta - z;
                // compute velocities
                const amrex::Real g = 9.81;
                const amrex::Real Omega =
                    std::sqrt(g * kappa * (1.0 + epsilon * epsilon));
                if (z < eta) {
                    vel(i, j, k, 0) = Omega * alpha * std::exp(kappa * z) *
                                      std::cos(kappa * x);
                    vel(i, j, k, 2) = Omega * alpha * std::exp(kappa * z) *
                                      std::sin(kappa * x);
                } else {
                    vel(i, j, k, 0) = vel_shear_rate * (z - eta);
                }
                // compute density
                amrex::Real smooth_heaviside;
                amrex::Real eps = std::cbrt(2. * dx[0] * dx[1] * dx[2]);
                if (phi(i, j, k) > eps) {
                    smooth_heaviside = 1.0;
                } else if (phi(i, j, k) < -eps) {
                    smooth_heaviside = 0.;
                } else {
                    smooth_heaviside =
                        0.5 * (1.0 + phi(i, j, k) / eps +
                               1.0 / utils::pi() *
                                   std::sin(phi(i, j, k) * utils::pi() / eps));
                }
                rho(i, j, k) =
                    rho1 * smooth_heaviside + rho2 * (1.0 - smooth_heaviside);
            });
    }
}

} // namespace amr_wind
