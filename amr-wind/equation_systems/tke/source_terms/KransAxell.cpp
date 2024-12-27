#include <AMReX_Orientation.H>

#include "amr-wind/equation_systems/tke/source_terms/KransAxell.H"
#include "amr-wind/CFDSim.H"
#include "amr-wind/turbulence/TurbulenceModel.H"
#include "amr-wind/utilities/linear_interpolation.H"
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
    , m_transport(sim.transport_model())
{
    AMREX_ALWAYS_ASSERT(sim.turbulence_model().model_name() == "KLAxell");
    auto coeffs = sim.turbulence_model().model_coeffs();
    amrex::ParmParse pp("ABL");
    pp.query("Cmu", m_Cmu);
    pp.query("kappa", m_kappa);
    pp.query("surface_roughness_z0", m_z0);
    pp.query("surface_temp_flux", m_heat_flux);
    pp.query("meso_sponge_start", m_sponge_start);
    {
        amrex::ParmParse pp_incflow("incflo");
        pp_incflow.queryarr("gravity", m_gravity);
    }
    m_ref_theta = m_transport.ref_theta();
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
    const auto& ref_theta_arr = (*m_ref_theta)(lev).const_array(mfi);
    const auto& geom = m_mesh.Geom(lev);
    const auto& problo = m_mesh.Geom(lev).ProbLoArray();
    const auto& probhi = m_mesh.Geom(lev).ProbHiArray();
    const auto& dx = geom.CellSizeArray();
    const auto& dt = m_time.delta_t();
    const amrex::Real heat_flux = m_heat_flux;
    const amrex::Real Cmu = m_Cmu;
    const amrex::Real sponge_start = m_sponge_start;
    const amrex::Real ref_tke = m_ref_tke;
    const auto tiny = std::numeric_limits<amrex::Real>::epsilon();
    const amrex::Real kappa = m_kappa;
    const amrex::Real z0 = m_z0;
    const bool has_terrain =
        this->m_sim.repo().int_field_exists("terrain_blank");
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> gravity{
        m_gravity[0], m_gravity[1], m_gravity[2]};
    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        amrex::Real bcforcing = 0;
        const amrex::Real ux = vel(i, j, k, 0);
        const amrex::Real uy = vel(i, j, k, 1);
        const amrex::Real z = problo[2] + (k + 0.5) * dx[2];
        if (k == 0) {
            const amrex::Real m = std::sqrt(ux * ux + uy * uy);
            const amrex::Real ustar = m * kappa / std::log(z / z0);
            const amrex::Real T0 = ref_theta_arr(i, j, k);
            const amrex::Real hf = std::abs(gravity[2]) / T0 * heat_flux;
            const amrex::Real rans_b = std::pow(
                std::max(hf, 0.0) * kappa * z / std::pow(Cmu, 3), (2.0 / 3.0));
            bcforcing =
                (ustar * ustar / (Cmu * Cmu) + rans_b - tke_arr(i, j, k)) / dt;
        }
        const amrex::Real zi =
            std::max((z - sponge_start) / (probhi[2] - sponge_start), 0.0);
        const amrex::Real sponge_forcing =
            zi * zi * (tke_arr(i, j, k) - ref_tke);
        dissip_arr(i, j, k) = std::pow(Cmu, 3) *
                              std::pow(tke_arr(i, j, k), 1.5) /
                              (tlscale_arr(i, j, k) + tiny);
        src_term(i, j, k) += shear_prod_arr(i, j, k) + buoy_prod_arr(i, j, k) -
                             dissip_arr(i, j, k) - sponge_forcing +
                             (1 - static_cast<int>(has_terrain)) * bcforcing;
    });
    // Add terrain components
    if (has_terrain) {
        const auto* const m_terrain_blank =
            &this->m_sim.repo().get_int_field("terrain_blank");
        const auto* const m_terrain_drag =
            &this->m_sim.repo().get_int_field("terrain_drag");
        const auto& blank_arr = (*m_terrain_blank)(lev).const_array(mfi);
        const auto& drag_arr = (*m_terrain_drag)(lev).const_array(mfi);
        auto* const m_terrain_height =
            &this->m_sim.repo().get_field("terrain_height");
        const auto& terrain_height = (*m_terrain_height)(lev).const_array(mfi);
        amrex::ParallelFor(
            bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                amrex::Real terrainforcing = 0;
                amrex::Real dragforcing = 0;
                const amrex::Real ux = vel(i, j, k, 0);
                const amrex::Real uy = vel(i, j, k, 1);
                amrex::Real z = 0.5 * dx[2];
                amrex::Real m = std::sqrt(ux * ux + uy * uy);
                const amrex::Real ustar = m * kappa / std::log(z / z0);
                const amrex::Real T0 = ref_theta_arr(i, j, k);
                const amrex::Real hf = std::abs(gravity[2]) / T0 * heat_flux;
                const amrex::Real rans_b = std::pow(
                    std::max(hf, 0.0) * kappa * z / std::pow(Cmu, 3),
                    (2.0 / 3.0));
                terrainforcing =
                    (ustar * ustar / (Cmu * Cmu) + rans_b - tke_arr(i, j, k)) /
                    dt;
                const amrex::Real uz = vel(i, j, k, 2);
                m = std::sqrt(ux * ux + uy * uy + uz * uz);
                const amrex::Real Cd =
                    std::min(10 / (dx[2] * m + tiny), 100 / dx[2]);
                dragforcing = -Cd * m * tke_arr(i, j, k, 0);
                z = std::max(
                    problo[2] + (k + 0.5) * dx[2] - terrain_height(i, j, k),
                    0.5 * dx[2]);
                const amrex::Real zi = std::max(
                    (z - sponge_start) / (probhi[2] - sponge_start), 0.0);
                const amrex::Real sponge_forcing =
                    zi * zi * (tke_arr(i, j, k) - ref_tke);
                src_term(i, j, k) += drag_arr(i, j, k) * terrainforcing +
                                     blank_arr(i, j, k) * dragforcing +
                                     static_cast<int>(has_terrain) * sponge_forcing;
            });
    }
}

} // namespace amr_wind::pde::tke
