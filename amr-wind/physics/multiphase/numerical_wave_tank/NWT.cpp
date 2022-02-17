#include "amr-wind/physics/multiphase/MultiPhase.H"
#include "amr-wind/physics/multiphase/numerical_wave_tank/NWT.H"
#include "amr-wind/physics/multiphase/numerical_wave_tank/wave_utils.H"
#include "amr-wind/physics/multiphase/numerical_wave_tank/wave_theories.H"
#include "amr-wind/utilities/trig_ops.H"
#include "amr-wind/CFDSim.H"
#include "AMReX_ParmParse.H"

namespace amr_wind {

NWT::NWT(CFDSim& sim)
    : m_sim(sim)
    , m_velocity(sim.repo().get_field("velocity"))
    , m_levelset(sim.repo().get_field("levelset"))
    , m_density(sim.repo().get_field("density"))
    , m_vof(sim.repo().get_field("vof"))
{
    amrex::ParmParse pp(identifier());
    pp.query("amplitude", m_amplitude);
    pp.query("wavelength", m_wavelength);
    pp.query("zero_sea_level", m_zsl);
    pp.query("water_depth", m_waterdepth);

    // Wave generation/absorption parameters
    pp.query("relax_zone_gen_length", m_gen_length);
    pp.query("numerical_beach_length", m_beach_length);
    pp.query("numerical_beach_start", m_x_start_beach);
}

NWT::~NWT() = default;

/** Initialize the velocity and vof fields at the beginning of the
 *  simulation.
 */
void NWT::initialize_fields(int level, const amrex::Geometry& geom)
{
    auto& velocity = m_velocity(level);
    auto& levelset = m_levelset(level);
    auto& density = m_density(level);
    velocity.setVal(0.0, 0, AMREX_SPACEDIM);

    const auto& dx = geom.CellSizeArray();
    const auto& problo = geom.ProbLoArray();

    const auto& mphase = m_sim.physics_manager().get<MultiPhase>();
    const amrex::Real rho1 = mphase.rho1();
    const amrex::Real rho2 = mphase.rho2();

    const amrex::Real zsl = m_zsl;

    for (amrex::MFIter mfi(levelset); mfi.isValid(); ++mfi) {
        const auto& vbx = mfi.growntilebox();
        auto vel = velocity.array(mfi);
        auto phi = levelset.array(mfi);
        auto rho = density.array(mfi);

        amrex::ParallelFor(
            vbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                const amrex::Real z = problo[2] + (k + 0.5) * dx[2];
                phi(i, j, k) = zsl - z;
                vel(i, j, k, 0) = 0.0;
                vel(i, j, k, 1) = 0.0;
                vel(i, j, k, 2) = 0.0;

                amrex::Real smooth_heaviside;
                amrex::Real eps = std::cbrt(2. * dx[0] * dx[1] * dx[2]);
                if (phi(i, j, k) > eps) {
                    smooth_heaviside = 1.0;
                } else if (phi(i, j, k) < -eps) {
                    smooth_heaviside = 0.;
                } else {
                    smooth_heaviside =
                        0.5 *
                        (1.0 + phi(i, j, k) / eps +
                         1.0 / M_PI * std::sin(phi(i, j, k) * M_PI / eps));
                }
                rho(i, j, k) =
                    rho1 * smooth_heaviside + rho2 * (1.0 - smooth_heaviside);
            });
    }
}

void NWT::pre_advance_work() { apply_relaxation_method(); }

void NWT::apply_relaxation_method()
{
    const auto& time = m_sim.time().current_time();
    const int nlevels = m_sim.repo().num_active_levels();
    const auto& geom = m_sim.mesh().Geom();

    amrex::Real ramp = (m_has_ramp) ? nwt::ramp(time, m_ramp_period) : 1.0;

    const auto& mphase = m_sim.physics_manager().get<MultiPhase>();
    const amrex::Real rho1 = mphase.rho1();
    const amrex::Real rho2 = mphase.rho2();

    // Interpolate within the relazation zones
    for (int lev = 0; lev < nlevels; ++lev) {
        for (amrex::MFIter mfi(m_vof(lev)); mfi.isValid(); ++mfi) {
            const auto& vbx = mfi.growntilebox();
            const auto& dx = geom[lev].CellSizeArray();
            const auto& problo = geom[lev].ProbLoArray();
            auto vel = m_velocity(lev).array(mfi);
            auto rho = m_density(lev).array(mfi);
            auto volfrac = m_vof(lev).array(mfi);
            const amrex::Real wavelength = m_wavelength;
            const amrex::Real waterdepth = m_waterdepth;
            const amrex::Real amplitude = m_amplitude;

            const amrex::Real gen_length = m_gen_length;
            const amrex::Real x_start_beach = m_x_start_beach;
            const amrex::Real beach_length = m_beach_length;
            const amrex::Real absorb_length_factor = m_absorb_length_factor;
            const amrex::Real zsl = m_zsl;

            amrex::ParallelFor(
                vbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    const amrex::Real x = problo[0] + (i + 0.5) * dx[0];
                    const amrex::Real y = problo[1] + (j + 0.5) * dx[1];
                    const amrex::Real z = problo[2] + (k + 0.5) * dx[2];
                    amrex::Vector<amrex::Real> wave_out(4);

                    wave_out = nwt::linear_monochromatic_waves(
                        wavelength, waterdepth, amplitude, x, y, z, time);

                    if (x <= gen_length) {
                        const amrex::Real Gamma =
                            nwt::Gamma_generate(x, gen_length);
                        const amrex::Real vf = (1 - Gamma) *
                                                   nwt::free_surface_to_vof(
                                                       wave_out[0], z, dx[2]) *
                                                   ramp +
                                               Gamma * volfrac(i, j, k);
                        // Doing clipping on the spot
                        volfrac(i, j, k) = (vf > 1. - 1.e-6) ? 1.0 : vf;
                        vel(i, j, k, 0) =
                            (1 - Gamma) * wave_out[1] * volfrac(i, j, k) *
                                ramp +
                            Gamma * vel(i, j, k, 0) * volfrac(i, j, k) +
                            (1. - volfrac(i, j, k)) * vel(i, j, k, 0);
                        vel(i, j, k, 1) =
                            (1 - Gamma) * wave_out[2] * volfrac(i, j, k) *
                                ramp +
                            Gamma * vel(i, j, k, 1) * volfrac(i, j, k) +
                            (1. - volfrac(i, j, k)) * vel(i, j, k, 1);
                        vel(i, j, k, 2) =
                            (1 - Gamma) * wave_out[3] * volfrac(i, j, k) *
                                ramp +
                            Gamma * vel(i, j, k, 2) * volfrac(i, j, k) +
                            (1. - volfrac(i, j, k)) * vel(i, j, k, 2);
                    }

                    if (x >= x_start_beach) {
                        const amrex::Real Gamma = nwt::Gamma_absorb(
                            x - x_start_beach, beach_length,
                            absorb_length_factor);
                        volfrac(i, j, k) =
                            (1.0 - Gamma) *
                                nwt::free_surface_to_vof(zsl, z, dx[2]) +
                            Gamma * volfrac(i, j, k);
                        vel(i, j, k, 0) =
                            Gamma * vel(i, j, k, 0) * volfrac(i, j, k);
                        vel(i, j, k, 1) =
                            Gamma * vel(i, j, k, 1) * volfrac(i, j, k);
                        vel(i, j, k, 2) =
                            Gamma * vel(i, j, k, 2) * volfrac(i, j, k);
                    }
                    // Make sure that density is updated before entering the
                    // solution
                    rho(i, j, k) = rho1 * volfrac(i, j, k) +
                                   rho2 * (1. - volfrac(i, j, k));
                });
        }
    }
} // namespace amr_wind

} // namespace amr_wind
