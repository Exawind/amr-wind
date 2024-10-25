#include <AMReX_Orientation.H>

#include "amr-wind/equation_systems/tke/source_terms/Krans.H"
#include "amr-wind/CFDSim.H"
#include "amr-wind/turbulence/TurbulenceModel.H"

namespace amr_wind::pde::tke {

Krans::Krans(const CFDSim& sim)
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
    AMREX_ALWAYS_ASSERT(sim.turbulence_model().model_name() == "OneEqRANS");
    auto coeffs = sim.turbulence_model().model_coeffs();
    amrex::ParmParse pp("ABL");
    pp.query("reference_temperature", m_ref_temp);
    pp.query("surface_temp_flux", m_heat_flux);
}

Krans::~Krans() = default;

void Krans::operator()(
    const int lev,
    const amrex::MFIter& mfi,
    const amrex::Box& bx,
    const FieldState fstate,
    const amrex::Array4<amrex::Real>& src_term) const
{
    const auto& vel =
        m_velocity.state(field_impl::dof_state(fstate))(lev).const_array(mfi);
    const bool has_terrain =
        this->m_sim.repo().int_field_exists("terrain_blank");
    const auto* const m_terrain_blank =
        has_terrain ? &this->m_sim.repo().get_int_field("terrain_blank")
                    : nullptr;
    const auto* const m_terrain_drag =
        has_terrain ? &this->m_sim.repo().get_int_field("terrain_drag")
                    : nullptr;
    const auto& blank_arr = has_terrain
                                ? (*m_terrain_blank)(lev).const_array(mfi)
                                : amrex::Array4<int>();
    const auto& drag_arr = has_terrain ? (*m_terrain_drag)(lev).const_array(mfi)
                                       : amrex::Array4<int>();
    const auto& tlscale_arr = (this->m_turb_lscale)(lev).array(mfi);
    const auto& shear_prod_arr = (this->m_shear_prod)(lev).array(mfi);
    const auto& buoy_prod_arr = (this->m_buoy_prod)(lev).array(mfi);
    const auto& dissip_arr = (this->m_dissip)(lev).array(mfi);
    const auto& tke_arr = m_tke(lev).array(mfi);
    const auto& geom = m_mesh.Geom(lev);
    const auto& problo = m_mesh.Geom(lev).ProbLoArray();
    const auto& dx = geom.CellSizeArray();
    const auto& dt = m_time.delta_t();
    const amrex::Real ref_temp = m_ref_temp;
    const amrex::Real heat_flux = 9.81 / ref_temp * m_heat_flux;
    const amrex::Real Cmu = 0.556;
    const auto tiny = std::numeric_limits<amrex::Real>::epsilon();
    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        amrex::Real bcforcing = 0;
        amrex::Real dragforcing = 0;
        amrex::Real ux = vel(i, j, k, 0);
        amrex::Real uy = vel(i, j, k, 1);
        const int cell_drag = (has_terrain) ? drag_arr(i, j, k) : 0;
        const int cell_blank = (has_terrain) ? blank_arr(i, j, k) : 0;
        amrex::Real ustar = 0.4;
        if (k == 0 && cell_drag == 1) {
            amrex::Real m = std::sqrt(ux * ux + uy * uy);
            const amrex::Real kappa = 0.41;
            const amrex::Real z0 = 0.1;
            const amrex::Real x3 =
                (k == 0) ? problo[2] + (k + 0.5) * dx[2] : 0.5 * dx[2];
            ustar = m * kappa / std::log(x3 / z0);
            const amrex::Real rans_b = std::pow(
                std::max(heat_flux, tiny) * 0.41 * x3 / std::pow(Cmu, 3),
                (2.0 / 3.0));
            bcforcing =
                (ustar * ustar / (0.556 * 0.556) + rans_b - tke_arr(i, j, k)) /
                dt;
        }
        if (cell_blank == 1) {
            amrex::Real uz = vel(i, j, k, 2);
            const amrex::Real m = std::sqrt(ux * ux + uy * uy + uz * uz);
            const amrex::Real Cd =
                std::min(10 / (dx[2] * m + tiny), 100 / dx[2]);
            dragforcing = -Cd * m * tke_arr(i, j, k, 0);
        }
        dissip_arr(i, j, k) = std::pow(Cmu, 3) *
                              std::pow(tke_arr(i, j, k), 1.5) /
                              (tlscale_arr(i, j, k) + tiny);
        src_term(i, j, k) += shear_prod_arr(i, j, k) + buoy_prod_arr(i, j, k) -
                             dissip_arr(i, j, k) + bcforcing + dragforcing;
    });
}

} // namespace amr_wind::pde::tke
