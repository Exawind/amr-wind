#include <AMReX_Orientation.H>

#include "amr-wind/equation_systems/tke/source_terms/KransAxell.H"
#include "amr-wind/CFDSim.H"
#include "amr-wind/turbulence/TurbulenceModel.H"

namespace amr_wind::pde::tke {

KransAxell::KransAxell(const CFDSim& sim)
    : m_turb_lscale(sim.repo().get_field("turb_lscale"))
    , m_shear_prod(sim.repo().get_field("shear_prod"))
    , m_buoy_prod(sim.repo().get_field("buoy_prod"))
    , m_dissip(sim.repo().get_field("dissipation"))
    , m_tke(sim.repo().get_field("tke"))
    , m_time(sim.time())
    , m_sim(sim)
    , m_mesh(sim.mesh())
    , m_velocity(sim.repo().get_field("velocity"))
{
    AMREX_ALWAYS_ASSERT(sim.turbulence_model().model_name() == "KLAxell");
    auto coeffs = sim.turbulence_model().model_coeffs();
    amrex::ParmParse pp("ABL");
    pp.query("reference_temperature", m_ref_temp);
    pp.query("surface_temp_flux", m_heat_flux);
    pp.query("meso_sponge_start", m_sponge_start);
}

KransAxell::~KransAxell() = default;

void KransAxell::operator()(
    const int lev,
    const amrex::MFIter& mfi,
    const amrex::Box& bx,
    const FieldState fstate,
    const amrex::Array4<amrex::Real>& src_term) const
{
    const auto& vel =
        m_velocity.state(field_impl::dof_state(fstate))(lev).const_array(mfi);
    const auto& tlscale_arr = (this->m_turb_lscale)(lev).array(mfi);
    const auto& shear_prod_arr = (this->m_shear_prod)(lev).array(mfi);
    const auto& buoy_prod_arr = (this->m_buoy_prod)(lev).array(mfi);
    const auto& dissip_arr = (this->m_dissip)(lev).array(mfi);
    const auto& tke_arr = m_tke(lev).array(mfi);
    const auto& geom = m_mesh.Geom(lev);
    const auto& problo = m_mesh.Geom(lev).ProbLoArray();
    const auto& probhi = m_mesh.Geom(lev).ProbHiArray();
    const auto& dx = geom.CellSizeArray();
    const auto& dt = m_time.delta_t();
    const amrex::Real ref_temp = m_ref_temp;
    const amrex::Real heat_flux = 9.81 / ref_temp * m_heat_flux;
    const amrex::Real Cmu = 0.556;
    const amrex::Real sponge_start = m_sponge_start;
    const amrex::Real ref_tke = m_ref_tke;
    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        amrex::Real bcforcing = 0;
        const amrex::Real ux = vel(i, j, k, 0);
        const amrex::Real uy = vel(i, j, k, 1);
        amrex::Real ustar = 0.4;
        if (k == 0) {
            const amrex::Real x3 = problo[2] + (k + 0.5) * dx[2];
            const amrex::Real m = std::sqrt(ux * ux + uy * uy);
            const amrex::Real kappa = 0.41;
            const amrex::Real z0 = 0.1;
            ustar = m * kappa / std::log(x3 / z0);
            const amrex::Real rans_b = std::pow(
                std::max(heat_flux, 0.0) * 0.41 * x3 / std::pow(0.556, 3),
                (2.0 / 3.0));
            bcforcing =
                (ustar * ustar / (0.556 * 0.556) + rans_b - tke_arr(i, j, k)) /
                dt;
        }
        const amrex::Real z = problo[2] + (k + 0.5) * dx[2];
        const amrex::Real zi =
            std::max((z - sponge_start) / (probhi[2] - sponge_start), 0.0);
        const amrex::Real sponge_forcing =
            zi * zi * (tke_arr(i, j, k) - ref_tke);
        dissip_arr(i, j, k) = std::pow(Cmu, 3) *
                              std::pow(tke_arr(i, j, k), 1.5) /
                              (tlscale_arr(i, j, k) + 1e-3);
        src_term(i, j, k) += shear_prod_arr(i, j, k) + buoy_prod_arr(i, j, k) -
                             dissip_arr(i, j, k) - sponge_forcing + bcforcing;
    });
    //! Add terrain components
    const bool has_terrain =
        this->m_sim.repo().int_field_exists("terrain_blank");
    if (has_terrain) {
        const auto* const m_terrain_blank =
            &this->m_sim.repo().get_int_field("terrain_blank");
        const auto* const m_terrain_drag =
            &this->m_sim.repo().get_int_field("terrain_drag");
        const auto& blank_arr = (*m_terrain_blank)(lev).const_array(mfi);
        const auto& drag_arr = (*m_terrain_drag)(lev).const_array(mfi);
        amrex::ParallelFor(
            bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                amrex::Real terrainforcing = 0;
                amrex::Real dragforcing = 0;
                const amrex::Real ux = vel(i, j, k, 0);
                const amrex::Real uy = vel(i, j, k, 1);
                amrex::Real ustar = 0.4;
                const amrex::Real x3 = 0.5 * dx[2];
                if (drag_arr(i, j, k) == 1) {
                    amrex::Real m = std::sqrt(ux * ux + uy * uy);
                    const amrex::Real kappa = 0.41;
                    const amrex::Real z0 = 0.1;
                    ustar = m * kappa / std::log(x3 / z0);
                    const amrex::Real rans_b = std::pow(
                        std::max(heat_flux, 0.0) * 0.41 * x3 /
                            std::pow(0.556, 3),
                        (2.0 / 3.0));
                    terrainforcing = (ustar * ustar / (0.556 * 0.556) + rans_b -
                                      tke_arr(i, j, k)) /
                                     dt;
                }
                if (blank_arr(i, j, k) == 1) {
                    const amrex::Real uz = vel(i, j, k, 2);
                    const amrex::Real m =
                        std::sqrt(ux * ux + uy * uy + uz * uz);
                    const amrex::Real Cd =
                        std::min(10 / (dx[2] * m + 1e-5), 100 / dx[2]);
                    dragforcing = -Cd * m * tke_arr(i, j, k, 0);
                }
                src_term(i, j, k) += terrainforcing + dragforcing;
            });
    }
}

} // namespace amr_wind::pde::tke
