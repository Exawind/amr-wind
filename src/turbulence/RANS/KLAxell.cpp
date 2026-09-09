#include "src/turbulence/RANS/KLAxell.H"
#include "src/equation_systems/PDEBase.H"
#include "src/turbulence/TurbModelDefs.H"
#include "src/fvm/gradient.H"
#include "src/fvm/strainrate.H"
#include "src/turbulence/turb_utils.H"
#include "src/equation_systems/tke/TKE.H"
#include "AMReX_ParmParse.H"
#include "src/utilities/math_ops.H"
#include "src/physics/ImmersedTerrain.H"
#include "src/physics/ImmersedWallModel.H"

using namespace amrex::literals;

namespace kynema_sgf {
namespace turbulence {

template <typename Transport>
KLAxell<Transport>::KLAxell(CFDSim& sim)
    : TurbModelBase<Transport>(sim)
    , m_vel(sim.repo().get_field("velocity"))
    , m_turb_lscale(sim.repo().declare_field("turb_lscale", 1))
    , m_shear_prod(sim.repo().declare_field("shear_prod", 1))
    , m_buoy_prod(sim.repo().declare_field("buoy_prod", 1))
    , m_dissip(sim.repo().declare_field("dissipation", 1))
    , m_rho(sim.repo().get_field("density"))
    , m_temperature(sim.repo().get_field("temperature"))
{
    auto& tke_eqn =
        sim.pde_manager().register_transport_pde(pde::TKE::pde_name());
    m_tke = &(tke_eqn.fields().field);
    auto& phy_mgr = this->m_sim.physics_manager();
    if (!phy_mgr.contains("ABL")) {
        amrex::Abort("KLAxell model only works with ABL physics");
    }
    {
        amrex::ParmParse pp("ABL");
        pp.get("surface_temp_flux", m_surf_flux);
        pp.query("meso_sponge_start", m_meso_sponge_start);
    }

    {
        amrex::ParmParse pp("incflo");
        pp.queryarr("gravity", m_gravity);
    }
    {
        amrex::ParmParse pp("KLAxell");
        pp.query("terrain_model", m_terrain_model);
        if (m_terrain_model != "TerrainDrag" &&
            m_terrain_model != "ImmersedTerrain") {
            amrex::Abort(
                "KLAxell.terrain_model must be TerrainDrag or ImmersedTerrain");
        }
        if (m_terrain_model == "ImmersedTerrain") {
            if (!sim.repo().field_exists("terrain_fraction")) {
                amrex::Abort(
                    "KLAxell.terrain_model = ImmersedTerrain requires the "
                    "ImmersedTerrain physics (terrain_fraction field not "
                    "found)");
            }
            amrex::ParmParse pp_terrain(
                immersedterrain::ImmersedTerrain::identifier());
            pp_terrain.query("solid_threshold", m_ib_solid_threshold);
            pp_terrain.query("drag_weight", m_ib_drag_weight);
        }
    }

    // TKE source term to be added to PDE
    turb_utils::inject_turbulence_src_terms(
        pde::TKE::pde_name(), {"KransAxell"});
}

template <typename Transport>
void KLAxell<Transport>::parse_model_coeffs()
{
    const std::string coeffs_dict = this->model_name() + "_coeffs";
    amrex::ParmParse pp(coeffs_dict);
    pp.query("Cmu", this->m_Cmu);
    pp.query("Cmu_prime", this->m_Cmu_prime);
    pp.query("Cb_stable", this->m_Cb_stable);
    pp.query("Cb_unstable", this->m_Cb_unstable);
    pp.query("prandtl", this->m_prandtl);
}

template <typename Transport>
TurbulenceModel::CoeffsDictType KLAxell<Transport>::model_coeffs() const
{
    return TurbulenceModel::CoeffsDictType{
        {"Cmu", this->m_Cmu},
        {"Cmu_prime", this->m_Cmu_prime},
        {"Cb_stable", this->m_Cb_stable},
        {"Cb_unstable", this->m_Cb_unstable},
        {"prandtl", this->m_prandtl}};
}

template <typename Transport>
void KLAxell<Transport>::post_init_actions()
{
    m_gradT = this->m_sim.repo().create_scratch_field(3, 0);
}

template <typename Transport>
void KLAxell<Transport>::post_regrid_actions()
{
    m_gradT = this->m_sim.repo().create_scratch_field(3, 0);
}

template <typename Transport>
void KLAxell<Transport>::update_turbulent_viscosity(
    const FieldState fstate, const DiffusionType /*unused*/)
{
    BL_PROFILE(
        "kynema-sgf::" + this->identifier() + "::update_turbulent_viscosity");

    fvm::gradient(*m_gradT, m_temperature.state(fstate));
    auto& gradT = *m_gradT;

    const auto& vel = this->m_vel.state(fstate);
    fvm::strainrate(this->m_shear_prod, vel);

    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> gravity{
        m_gravity[0], m_gravity[1], m_gravity[2]};
    const auto beta = (this->m_transport).beta();
    const amrex::Real Cmu = m_Cmu;
    const amrex::Real Cb_stable = m_Cb_stable;
    const amrex::Real Cb_unstable = m_Cb_unstable;
    auto& mu_turb = this->mu_turb();
    const auto& den = this->m_rho.state(fstate);
    const auto& repo = mu_turb.repo();
    const auto& geom_vec = repo.mesh().Geom();
    const int nlevels = repo.num_active_levels();
    const amrex::Real Rtc = -1.0_rt;
    const amrex::Real Rtmin = -3.0_rt;
    const amrex::Real lambda = 30.0_rt;
    const amrex::Real kappa = 0.41_rt;
    const amrex::Real surf_flux = m_surf_flux;
    const auto tiny = std::numeric_limits<amrex::Real>::epsilon();
    const amrex::Real lengthscale_switch = m_meso_sponge_start;
    const bool use_immersed = (m_terrain_model == "ImmersedTerrain");
    for (int lev = 0; lev < nlevels; ++lev) {
        if (use_immersed) {
            immersed_terrain_viscosity(lev, fstate, gradT, *beta);
            continue;
        }
        const auto& geom = geom_vec[lev];
        const auto& problo = repo.mesh().Geom(lev).ProbLoArray();
        const amrex::Real dz = geom.CellSize()[2];

        const auto& mu_arrs = mu_turb(lev).arrays();
        const auto& rho_arrs = den(lev).const_arrays();
        const auto& gradT_arrs = gradT(lev).const_arrays();
        const auto& tlscale_arrs = (this->m_turb_lscale)(lev).arrays();
        const auto& tke_arrs = (*this->m_tke)(lev).arrays();
        const auto& buoy_prod_arrs = (this->m_buoy_prod)(lev).arrays();
        const auto& shear_prod_arrs = (this->m_shear_prod)(lev).arrays();
        const auto& beta_arrs = (*beta)(lev).const_arrays();

        //! Add terrain components
        const bool has_terrain =
            this->m_sim.repo().int_field_exists("terrain_blank");
        if (has_terrain) {
            const auto* m_terrain_height =
                &this->m_sim.repo().get_field("terrain_height");
            const auto* m_terrain_blank =
                &this->m_sim.repo().get_int_field("terrain_blank");
            const auto& ht_arrs = (*m_terrain_height)(lev).const_arrays();
            const auto& blank_arrs = (*m_terrain_blank)(lev).const_arrays();
            amrex::ParallelFor(
                mu_turb(lev),
                [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
                    amrex::Real stratification =
                        -((gradT_arrs[nbx](i, j, k, 0) * gravity[0]) +
                          (gradT_arrs[nbx](i, j, k, 1) * gravity[1]) +
                          (gradT_arrs[nbx](i, j, k, 2) * gravity[2])) *
                        beta_arrs[nbx](i, j, k);
                    const amrex::Real z = amrex::max<amrex::Real>(
                        problo[2] + ((k + 0.5_rt) * dz) - ht_arrs[nbx](i, j, k),
                        0.5_rt * dz);
                    const amrex::Real lscale_s =
                        (lambda * kappa * z) / (lambda + kappa * z);
                    const amrex::Real lscale_b =
                        Cb_stable *
                        std::sqrt(
                            tke_arrs[nbx](i, j, k) /
                            amrex::max<amrex::Real>(stratification, tiny));
                    amrex::Real epsilon =
                        utils::powi(Cmu, 3) *
                        std::pow(tke_arrs[nbx](i, j, k), 1.5_rt) /
                        (tlscale_arrs[nbx](i, j, k) + tiny);
                    amrex::Real Rt =
                        utils::powi(tke_arrs[nbx](i, j, k) / epsilon, 2) *
                        stratification;
                    Rt = (Rt > Rtc)
                             ? Rt
                             : amrex::max<amrex::Real>(
                                   Rt, Rt - (utils::powi(Rt - Rtc, 2) /
                                             (Rt + Rtmin - 2.0_rt * Rtc)));
                    tlscale_arrs[nbx](i, j, k) =
                        (stratification > 0)
                            ? std::sqrt(
                                  utils::powi(lscale_s * lscale_b, 2) /
                                  (utils::powi(lscale_s, 2) +
                                   utils::powi(lscale_b, 2)))
                            : lscale_s *
                                  std::sqrt(
                                      1.0_rt -
                                      (utils::powi(Cmu, 6) *
                                       utils::powi(Cb_unstable, -2) * Rt));
                    tlscale_arrs[nbx](i, j, k) =
                        (stratification > 0)
                            ? amrex::min<amrex::Real>(
                                  tlscale_arrs[nbx](i, j, k),
                                  std::sqrt(
                                      Cmu * tke_arrs[nbx](i, j, k) /
                                      stratification))
                            : tlscale_arrs[nbx](i, j, k);
                    tlscale_arrs[nbx](i, j, k) =
                        (std::abs(surf_flux) < 1.0e-5_rt &&
                         z <= lengthscale_switch)
                            ? lscale_s
                            : tlscale_arrs[nbx](i, j, k);
                    Rt = (std::abs(surf_flux) < 1.0e-5_rt &&
                          z <= lengthscale_switch)
                             ? 0.0_rt
                             : Rt;
                    const amrex::Real Cmu_Rt =
                        (Cmu + 0.108_rt * Rt) /
                        (1.0_rt + 0.308_rt * Rt +
                         0.00837_rt * utils::powi(Rt, 2));
                    mu_arrs[nbx](i, j, k) = rho_arrs[nbx](i, j, k) * Cmu_Rt *
                                            tlscale_arrs[nbx](i, j, k) *
                                            std::sqrt(tke_arrs[nbx](i, j, k)) *
                                            (1.0_rt - blank_arrs[nbx](i, j, k));
                    const amrex::Real Cmu_prime_Rt =
                        Cmu / (1.0_rt + 0.277_rt * Rt);
                    const amrex::Real muPrime =
                        rho_arrs[nbx](i, j, k) * Cmu_prime_Rt *
                        tlscale_arrs[nbx](i, j, k) *
                        std::sqrt(tke_arrs[nbx](i, j, k)) *
                        (1.0_rt - blank_arrs[nbx](i, j, k));
                    buoy_prod_arrs[nbx](i, j, k) = -muPrime * stratification;
                    shear_prod_arrs[nbx](i, j, k) *=
                        shear_prod_arrs[nbx](i, j, k) * mu_arrs[nbx](i, j, k);
                });
        } else {
            amrex::ParallelFor(
                mu_turb(lev),
                [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
                    amrex::Real stratification =
                        -((gradT_arrs[nbx](i, j, k, 0) * gravity[0]) +
                          (gradT_arrs[nbx](i, j, k, 1) * gravity[1]) +
                          (gradT_arrs[nbx](i, j, k, 2) * gravity[2])) *
                        beta_arrs[nbx](i, j, k);
                    const amrex::Real z = problo[2] + ((k + 0.5_rt) * dz);
                    const amrex::Real lscale_s =
                        (lambda * kappa * z) / (lambda + kappa * z);
                    const amrex::Real lscale_b =
                        Cb_stable *
                        std::sqrt(
                            tke_arrs[nbx](i, j, k) /
                            amrex::max<amrex::Real>(stratification, tiny));
                    amrex::Real epsilon =
                        utils::powi(Cmu, 3) *
                        std::pow(tke_arrs[nbx](i, j, k), 1.5_rt) /
                        (tlscale_arrs[nbx](i, j, k) + tiny);
                    amrex::Real Rt =
                        utils::powi(tke_arrs[nbx](i, j, k) / epsilon, 2) *
                        stratification;
                    Rt = (Rt > Rtc)
                             ? Rt
                             : amrex::max<amrex::Real>(
                                   Rt, Rt - (utils::powi(Rt - Rtc, 2) /
                                             (Rt + Rtmin - 2.0_rt * Rtc)));
                    tlscale_arrs[nbx](i, j, k) =
                        (stratification > 0)
                            ? std::sqrt(
                                  utils::powi(lscale_s * lscale_b, 2) /
                                  (utils::powi(lscale_s, 2) +
                                   utils::powi(lscale_b, 2)))
                            : lscale_s *
                                  std::sqrt(
                                      1.0_rt -
                                      (utils::powi(Cmu, 6) *
                                       utils::powi(Cb_unstable, -2) * Rt));
                    tlscale_arrs[nbx](i, j, k) =
                        (stratification > 0)
                            ? amrex::min<amrex::Real>(
                                  tlscale_arrs[nbx](i, j, k),
                                  std::sqrt(
                                      Cmu * tke_arrs[nbx](i, j, k) /
                                      stratification))
                            : tlscale_arrs[nbx](i, j, k);
                    tlscale_arrs[nbx](i, j, k) =
                        (std::abs(surf_flux) < 1.0e-5_rt &&
                         z <= lengthscale_switch)
                            ? lscale_s
                            : tlscale_arrs[nbx](i, j, k);
                    Rt = (std::abs(surf_flux) < 1.0e-5_rt &&
                          z <= lengthscale_switch)
                             ? 0.0_rt
                             : Rt;
                    const amrex::Real Cmu_Rt =
                        (Cmu + 0.108_rt * Rt) /
                        (1.0_rt + 0.308_rt * Rt +
                         0.00837_rt * utils::powi(Rt, 2));
                    mu_arrs[nbx](i, j, k) = rho_arrs[nbx](i, j, k) * Cmu_Rt *
                                            tlscale_arrs[nbx](i, j, k) *
                                            std::sqrt(tke_arrs[nbx](i, j, k));
                    const amrex::Real Cmu_prime_Rt =
                        Cmu / (1.0_rt + 0.277_rt * Rt);
                    const amrex::Real muPrime =
                        rho_arrs[nbx](i, j, k) * Cmu_prime_Rt *
                        tlscale_arrs[nbx](i, j, k) *
                        std::sqrt(tke_arrs[nbx](i, j, k));
                    buoy_prod_arrs[nbx](i, j, k) = -muPrime * stratification;
                    shear_prod_arrs[nbx](i, j, k) *=
                        shear_prod_arrs[nbx](i, j, k) * mu_arrs[nbx](i, j, k);
                });
        }
    }
    amrex::Gpu::streamSynchronize();

    mu_turb.fillpatch(this->m_sim.time().current_time());
}

// KLAxell viscosity with the ImmersedTerrain fields: the TerrainDrag kernel
//  with the height above the terrain taken from ``terrain_surface`` and the
//  binary blanking replaced by (1 - w_solid), w_solid being the drag weight
//  of the cell (``ImmersedTerrain.drag_weight``), so that partial cells keep
//  a fraction of the eddy viscosity.
template <typename Transport>
void KLAxell<Transport>::immersed_terrain_viscosity(
    const int lev,
    const FieldState fstate,
    const ScratchField& gradT,
    const ScratchField& beta_field)
{
    using immersedterrain::ImmersedTerrain;

    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> gravity{
        m_gravity[0], m_gravity[1], m_gravity[2]};
    const amrex::Real Cmu = m_Cmu;
    const amrex::Real Cb_stable = m_Cb_stable;
    const amrex::Real Cb_unstable = m_Cb_unstable;
    const amrex::Real Rtc = -1.0_rt;
    const amrex::Real Rtmin = -3.0_rt;
    const amrex::Real lambda = 30.0_rt;
    const amrex::Real kappa = 0.41_rt;
    const amrex::Real surf_flux = m_surf_flux;
    const auto tiny = std::numeric_limits<amrex::Real>::epsilon();
    const amrex::Real lengthscale_switch = m_meso_sponge_start;
    const bool center_weight = (m_ib_drag_weight == "center");
    const amrex::Real solid_threshold = m_ib_solid_threshold;

    auto& mu_turb = this->mu_turb();
    const auto& den = this->m_rho.state(fstate);
    const auto& repo = mu_turb.repo();
    const auto& geom = repo.mesh().Geom(lev);
    const auto& problo = geom.ProbLoArray();
    const amrex::Real dz = geom.CellSize()[2];

    const auto& mu_arrs = mu_turb(lev).arrays();
    const auto& rho_arrs = den(lev).const_arrays();
    const auto& gradT_arrs = gradT(lev).const_arrays();
    const auto& tlscale_arrs = (this->m_turb_lscale)(lev).arrays();
    const auto& tke_arrs = (*this->m_tke)(lev).arrays();
    const auto& buoy_prod_arrs = (this->m_buoy_prod)(lev).arrays();
    const auto& shear_prod_arrs = (this->m_shear_prod)(lev).arrays();
    const auto& beta_arrs = beta_field(lev).const_arrays();
    const auto& frac_arrs =
        this->m_sim.repo().get_field("terrain_fraction")(lev).const_arrays();
    const auto& surf_arrs =
        this->m_sim.repo().get_field("terrain_surface")(lev).const_arrays();

    amrex::ParallelFor(
        mu_turb(lev), [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
            const amrex::Real w_solid = immersed_wall::solid_weight(
                frac_arrs[nbx](i, j, k, 0), center_weight, solid_threshold);
            const amrex::Real fluid_weight = 1.0_rt - w_solid;
            const amrex::Real h =
                surf_arrs[nbx](i, j, k, ImmersedTerrain::surf_height);
            amrex::Real stratification =
                -((gradT_arrs[nbx](i, j, k, 0) * gravity[0]) +
                  (gradT_arrs[nbx](i, j, k, 1) * gravity[1]) +
                  (gradT_arrs[nbx](i, j, k, 2) * gravity[2])) *
                beta_arrs[nbx](i, j, k);
            const amrex::Real z = amrex::max<amrex::Real>(
                problo[2] + ((k + 0.5_rt) * dz) - h, 0.5_rt * dz);
            const amrex::Real lscale_s =
                (lambda * kappa * z) / (lambda + (kappa * z));
            const amrex::Real lscale_b =
                Cb_stable * std::sqrt(
                                tke_arrs[nbx](i, j, k) /
                                amrex::max<amrex::Real>(stratification, tiny));
            const amrex::Real epsilon =
                utils::powi(Cmu, 3) * std::pow(tke_arrs[nbx](i, j, k), 1.5_rt) /
                (tlscale_arrs[nbx](i, j, k) + tiny);
            amrex::Real Rt = utils::powi(tke_arrs[nbx](i, j, k) / epsilon, 2) *
                             stratification;
            Rt = (Rt > Rtc) ? Rt
                            : amrex::max<amrex::Real>(
                                  Rt, Rt - (utils::powi(Rt - Rtc, 2) /
                                            (Rt + Rtmin - (2.0_rt * Rtc))));
            tlscale_arrs[nbx](i, j, k) =
                (stratification > 0)
                    ? std::sqrt(
                          utils::powi(lscale_s * lscale_b, 2) /
                          (utils::powi(lscale_s, 2) + utils::powi(lscale_b, 2)))
                    : lscale_s *
                          std::sqrt(
                              1.0_rt - (utils::powi(Cmu, 6) *
                                        utils::powi(Cb_unstable, -2) * Rt));
            tlscale_arrs[nbx](i, j, k) =
                (stratification > 0)
                    ? amrex::min<amrex::Real>(
                          tlscale_arrs[nbx](i, j, k),
                          std::sqrt(
                              Cmu * tke_arrs[nbx](i, j, k) / stratification))
                    : tlscale_arrs[nbx](i, j, k);
            const bool neutral_switch =
                (std::abs(surf_flux) < 1.0e-5_rt) && (z <= lengthscale_switch);
            tlscale_arrs[nbx](i, j, k) =
                neutral_switch ? lscale_s : tlscale_arrs[nbx](i, j, k);
            Rt = neutral_switch ? 0.0_rt : Rt;
            const amrex::Real Cmu_Rt =
                (Cmu + (0.108_rt * Rt)) /
                (1.0_rt + (0.308_rt * Rt) + (0.00837_rt * utils::powi(Rt, 2)));
            mu_arrs[nbx](i, j, k) =
                rho_arrs[nbx](i, j, k) * Cmu_Rt * tlscale_arrs[nbx](i, j, k) *
                std::sqrt(tke_arrs[nbx](i, j, k)) * fluid_weight;
            const amrex::Real Cmu_prime_Rt = Cmu / (1.0_rt + (0.277_rt * Rt));
            const amrex::Real muPrime = rho_arrs[nbx](i, j, k) * Cmu_prime_Rt *
                                        tlscale_arrs[nbx](i, j, k) *
                                        std::sqrt(tke_arrs[nbx](i, j, k)) *
                                        fluid_weight;
            buoy_prod_arrs[nbx](i, j, k) = -muPrime * stratification;
            shear_prod_arrs[nbx](i, j, k) *=
                shear_prod_arrs[nbx](i, j, k) * mu_arrs[nbx](i, j, k);
        });
}

template <typename Transport>
void KLAxell<Transport>::update_alphaeff(Field& alphaeff)
{

    BL_PROFILE("kynema-sgf::" + this->identifier() + "::update_alphaeff");
    auto lam_alpha = (this->m_transport).alpha();
    auto& mu_turb = this->m_mu_turb;
    auto& repo = mu_turb.repo();

    fvm::gradient(*m_gradT, m_temperature);
    auto& gradT = *m_gradT;
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> gravity{
        m_gravity[0], m_gravity[1], m_gravity[2]};
    const auto beta = (this->m_transport).beta();
    const amrex::Real Cmu = m_Cmu;
    const int nlevels = repo.num_active_levels();
    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& muturb_arrs = mu_turb(lev).arrays();
        const auto& alphaeff_arrs = alphaeff(lev).arrays();
        const auto& lam_diff_arrs = (*lam_alpha)(lev).arrays();
        const auto& tke_arrs = (*this->m_tke)(lev).arrays();
        const auto& gradT_arrs = gradT(lev).const_arrays();
        const auto& tlscale_arrs = (this->m_turb_lscale)(lev).arrays();
        const auto& beta_arrs = (*beta)(lev).const_arrays();
        const amrex::Real Rtc = -1.0_rt;
        const amrex::Real Rtmin = -3.0_rt;
        amrex::ParallelFor(
            mu_turb(lev), [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
                amrex::Real stratification =
                    -((gradT_arrs[nbx](i, j, k, 0) * gravity[0]) +
                      (gradT_arrs[nbx](i, j, k, 1) * gravity[1]) +
                      (gradT_arrs[nbx](i, j, k, 2) * gravity[2])) *
                    beta_arrs[nbx](i, j, k);
                amrex::Real epsilon = utils::powi(Cmu, 3) *
                                      std::pow(tke_arrs[nbx](i, j, k), 1.5_rt) /
                                      tlscale_arrs[nbx](i, j, k);
                amrex::Real Rt =
                    utils::powi(tke_arrs[nbx](i, j, k) / epsilon, 2) *
                    stratification;
                Rt = (Rt > Rtc) ? Rt
                                : amrex::max<amrex::Real>(
                                      Rt, Rt - (utils::powi(Rt - Rtc, 2) /
                                                (Rt + Rtmin - 2.0_rt * Rtc)));
                const amrex::Real prandtlRt =
                    (1.0_rt + 0.193_rt * Rt) / (1.0_rt + 0.0302_rt * Rt);
                alphaeff_arrs[nbx](i, j, k) =
                    lam_diff_arrs[nbx](i, j, k) +
                    (muturb_arrs[nbx](i, j, k) / prandtlRt);
            });
    }
    amrex::Gpu::streamSynchronize();

    alphaeff.fillpatch(this->m_sim.time().current_time());
}

template <typename Transport>
void KLAxell<Transport>::update_scalar_diff(
    Field& deff, const std::string& name)
{
    BL_PROFILE("kynema-sgf::" + this->identifier() + "::update_scalar_diff");

    if (name == pde::TKE::var_name()) {
        auto& mu_turb = this->mu_turb();
        deff.setVal(0.0_rt);
        field_ops::saxpy(
            deff, 2.0_rt, mu_turb, 0, 0, deff.num_comp(), deff.num_grow());
    } else {
        amrex::Abort(
            "KLAxell:update_scalar_diff not implemented for field " + name);
    }
}

template <typename Transport>
void KLAxell<Transport>::post_advance_work()
{
    BL_PROFILE("kynema-sgf::" + this->identifier() + "::post_advance_work");
}

} // namespace turbulence

INSTANTIATE_TURBULENCE_MODEL(KLAxell);

} // namespace kynema_sgf
