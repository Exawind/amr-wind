#include <cmath>
#include "src/turbulence/LES/Kosovic.H"
#include "src/turbulence/TurbModelDefs.H"
#include "src/fvm/nonLinearSum.H"
#include "src/fvm/strainrate.H"
#include "src/fvm/divergence.H"
#include "src/utilities/math_ops.H"
#include "src/fvm/gradient.H"
#include "AMReX_REAL.H"
#include "AMReX_MultiFab.H"
#include "AMReX_ParmParse.H"
#include "src/wind_energy/ABL.H"
#include "src/physics/ImmersedTerrain.H"
#include "src/physics/ImmersedWallModel.H"

using namespace amrex::literals;

namespace kynema_sgf {
namespace turbulence {

template <typename Transport>
// cppcheck-suppress uninitMemberVar
Kosovic<Transport>::Kosovic(CFDSim& sim)
    : TurbModelBase<Transport>(sim)
    , m_vel(sim.repo().get_field("velocity"))
    , m_rho(sim.repo().get_field("density"))
    , m_theta(sim.repo().get_field("temperature"))
    , m_Nij(sim.repo().declare_field("Nij", 9, 1, 1))
    , m_divNij(sim.repo().declare_field("divNij", 3))
{
    amrex::ParmParse pp("Kosovic");
    pp.query("Cb", m_Cb);
    m_Cs = std::sqrt(
        8.0_rt * (1.0_rt + m_Cb) /
        (27.0_rt * std::numbers::pi_v<amrex::Real> *
         std::numbers::pi_v<amrex::Real>));
    m_C1 = std::sqrt(960.0_rt) * m_Cb / (7.0_rt * (1.0_rt + m_Cb) * m_Sk);
    m_C2 = m_C1;
    pp.query("surfaceRANS", m_surfaceRANS);
    if (m_surfaceRANS) {
        m_surfaceFactor = 1.0_rt;
        pp.query("switchLoc", m_switchLoc);
        pp.query("surfaceRANSExp", m_surfaceRANSExp);
    } else {
        m_surfaceFactor = 0.0_rt;
    }
    pp.query("writeTerms", m_writeTerms);
    if (m_writeTerms) {
        this->m_sim.io_manager().register_io_var("Nij");
        this->m_sim.io_manager().register_io_var("divNij");
    }
    pp.query("LESOff", m_LESTurnOff);
    pp.query("muCoeff", m_muCoeff);
    amrex::ParmParse pp_abl("ABL");
    pp_abl.query("wall_het_model", m_wall_het_model);
    pp_abl.query("monin_obukhov_length", m_monin_obukhov_length);
    pp_abl.query("kappa", m_kappa);
    pp_abl.query("mo_gamma_m", m_gamma_m);
    pp_abl.query("mo_beta_m", m_beta_m);
    pp_abl.query("surface_roughness_z0", m_surface_roughness_z0);

    amrex::ParmParse pp_incflo("incflo");
    pp_incflo.queryarr("gravity", m_gravity);

    pp.query("terrain_model", m_terrain_model);
    if (m_terrain_model != "TerrainDrag" &&
        m_terrain_model != "ImmersedTerrain") {
        amrex::Abort(
            "Kosovic.terrain_model must be TerrainDrag or ImmersedTerrain");
    }
    if (m_terrain_model == "ImmersedTerrain") {
        if (!sim.repo().field_exists("terrain_fraction")) {
            amrex::Abort(
                "Kosovic.terrain_model = ImmersedTerrain requires the "
                "ImmersedTerrain physics (terrain_fraction field not found)");
        }
        // Same wall-model and surface-cell settings as the source terms so
        // that the SGS viscosity and the wall forcing see one friction
        // velocity
        amrex::ParmParse pp_ib("ImmersedDragForcing");
        pp_ib.query("wall_model", m_ib_wall_model);
        pp_ib.query("reference_distance", m_ib_reference_distance);
        pp_ib.query("minimum_z0", m_ib_min_z0);
        amrex::ParmParse pp_terrain(
            immersedterrain::ImmersedTerrain::identifier());
        pp_terrain.query("solid_threshold", m_ib_solid_threshold);
        pp_terrain.query("drag_weight", m_ib_drag_weight);
    }
}
template <typename Transport>
void Kosovic<Transport>::update_turbulent_viscosity(
    const FieldState fstate, const DiffusionType /*unused*/)
{
    BL_PROFILE(
        "kynema-sgf::" + this->identifier() + "::update_turbulent_viscosity");

    auto& mu_turb = this->mu_turb();
    const auto& repo = mu_turb.repo();
    const auto& vel = m_vel.state(fstate);
    const auto& den = m_rho.state(fstate);
    const auto ref_theta = (this->m_transport).ref_theta();

    auto gradT = (this->m_sim.repo()).create_scratch_field(3, 0);
    fvm::gradient(*gradT, m_theta.state(fstate));
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> gravity{
        m_gravity[0], m_gravity[1], m_gravity[2]};

    const amrex::Real Cs_sqr = this->m_Cs * this->m_Cs;
    const amrex::Real monin_obukhov_length = m_monin_obukhov_length;
    const amrex::Real kappa = m_kappa;
    const amrex::Real surface_roughness_z0 = m_surface_roughness_z0;
    const amrex::Real z0_min = 1.0e-4_rt;
    const amrex::Real dMdz_min = 0.01_rt;
    const amrex::Real locLESTurnOff = m_LESTurnOff;
    const amrex::Real locSwitchLoc = m_switchLoc;
    const amrex::Real locSurfaceRANSExp = m_surfaceRANSExp;
    const amrex::Real locSurfaceFactor = m_surfaceFactor;
    const amrex::Real locC1 = m_C1;
    const amrex::Real tol = constants::TIGHT_TOL;

    const auto& geom_vec = repo.mesh().Geom();
    // TerrainDrag fields (binary blanking, drag cells); the ImmersedTerrain
    // fields are handled in a separate kernel selected by terrain_model
    const bool use_immersed = (m_terrain_model == "ImmersedTerrain");
    const bool has_terrain =
        !use_immersed && this->m_sim.repo().int_field_exists("terrain_blank");
    const auto* m_terrain_blank =
        has_terrain ? &this->m_sim.repo().get_int_field("terrain_blank")
                    : nullptr;
    const auto* m_terrain_drag =
        has_terrain ? &this->m_sim.repo().get_int_field("terrain_drag")
                    : nullptr;
    const auto* m_terrain_height =
        has_terrain ? &this->m_sim.repo().get_field("terrain_height") : nullptr;
    const auto* m_terrain_z0 =
        has_terrain ? &this->m_sim.repo().get_field("terrainz0") : nullptr;
    // Populate strainrate into the turbulent viscosity arrays to avoid creating
    // a temporary buffer
    fvm::strainrate(mu_turb, vel);
    // Non-linear component Nij is computed here and goes into Body Forcing
    fvm::nonlinearsum(m_Nij, vel);
    fvm::divergence(m_divNij, m_Nij);
    const int nlevels = repo.num_active_levels();
    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& geom = geom_vec[lev];
        const auto& problo = repo.mesh().Geom(lev).ProbLoArray();
        const amrex::Real dx = geom.CellSize()[0];
        const amrex::Real dy = geom.CellSize()[1];
        const amrex::Real dz = geom.CellSize()[2];
        const amrex::Real ds = std::cbrt(dx * dy * dz);
        const amrex::Real ds_sqr = ds * ds;
        const amrex::Real smag_factor = Cs_sqr * ds_sqr;
        const auto& mu_arrs = mu_turb(lev).arrays();
        const auto& rho_arrs = den(lev).const_arrays();
        const auto& vel_arrs = vel(lev).const_arrays();
        const auto& gradT_arrs = (*gradT)(lev).const_arrays();
        const auto& divNij_arrs = (this->m_divNij)(lev).arrays();
        const auto& blank_arrs = has_terrain
                                     ? (*m_terrain_blank)(lev).const_arrays()
                                     : amrex::MultiArray4<const int>();
        const auto& drag_arrs = has_terrain
                                    ? (*m_terrain_drag)(lev).const_arrays()
                                    : amrex::MultiArray4<const int>();
        const auto& height_arrs = has_terrain
                                      ? (*m_terrain_height)(lev).const_arrays()
                                      : amrex::MultiArray4<const amrex::Real>();
        const auto& z0_arrs = has_terrain
                                  ? (*m_terrain_z0)(lev).const_arrays()
                                  : amrex::MultiArray4<const amrex::Real>();
        const amrex::Real non_neutral_neighbour =
            (m_wall_het_model == "mol")
                ? MOData::calc_psi_m(
                      1.5_rt * dz / monin_obukhov_length, m_beta_m, m_gamma_m)
                : 0.0_rt;
        const auto& ref_theta_arrs = (*ref_theta)(lev).const_arrays();
        if (use_immersed) {
            immersed_terrain_viscosity(lev, fstate, *gradT, *ref_theta);
            continue;
        }
        amrex::ParallelFor(
            mu_turb(lev), [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
                const amrex::Real rho = rho_arrs[nbx](i, j, k);
                amrex::Real x3 = problo[2] + ((k + 0.5_rt) * dz);
                x3 = (has_terrain)
                         ? amrex::max<amrex::Real>(
                               x3 - height_arrs[nbx](i, j, k, 0), 0.5_rt * dz)
                         : x3;
                const amrex::Real fmu = std::exp(-x3 / locSwitchLoc);
                const amrex::Real phiM =
                    (monin_obukhov_length < 0)
                        ? std::pow(
                              1.0_rt - (16.0_rt * x3 / monin_obukhov_length),
                              -0.25_rt)
                        : 1.0_rt + (5.0_rt * x3 / monin_obukhov_length);
                const amrex::Real wall_distance =
                    (has_terrain)
                        ? amrex::max<amrex::Real>(
                              ((k + 1) * dz) - height_arrs[nbx](i, j, k, 0), dz)
                        : (k + 1) * dz;
                const amrex::Real ransL =
                    utils::powi(0.41_rt * wall_distance / phiM, 2);
                const amrex::Real turnOff = std::exp(-x3 / locLESTurnOff);
                const amrex::Real viscosityScale =
                    (locSurfaceFactor *
                     (std::pow(1.0_rt - fmu, locSurfaceRANSExp) * smag_factor +
                      std::pow(fmu, locSurfaceRANSExp) * ransL)) +
                    ((1.0_rt - locSurfaceFactor) * smag_factor);
                const amrex::Real blankTerrain =
                    (has_terrain) ? 1 - blank_arrs[nbx](i, j, k, 0) : 1.0_rt;
                const amrex::Real mut =
                    mu_arrs[nbx](i, j, k) * mu_arrs[nbx](i, j, k);
                const amrex::Real T0 = ref_theta_arrs[nbx](i, j, k);
                const amrex::Real stratification_sensor =
                    -(gradT_arrs[nbx](i, j, k, 0) * gravity[0] +
                      gradT_arrs[nbx](i, j, k, 1) * gravity[1] +
                      gradT_arrs[nbx](i, j, k, 2) * gravity[2]) /
                    T0;
                amrex::Real stratification = 1.0_rt;
                amrex::Real non_linear_coeff = 1.0_rt;

                if (stratification_sensor > tol) {
                    // stable
                    non_linear_coeff =
                        (mut - 3.0_rt * stratification_sensor < tol) ? 0.0_rt
                                                                     : 1.0_rt;
                    stratification = std::sqrt(
                        amrex::max<amrex::Real>(
                            tol, mut - 3.0_rt * stratification_sensor));
                } else {
                    stratification = std::sqrt(mut);
                }

                mu_arrs[nbx](i, j, k) = rho * viscosityScale * turnOff *
                                        blankTerrain * stratification;
                // log-law
                const amrex::Real ux = vel_arrs[nbx](i, j, k + 1, 0);
                const amrex::Real uy = vel_arrs[nbx](i, j, k + 1, 1);
                const amrex::Real m = std::sqrt((ux * ux) + (uy * uy));
                const amrex::Real local_z0 =
                    (has_terrain) ? amrex::max<amrex::Real>(
                                        z0_arrs[nbx](i, j, k, 0), z0_min)
                                  : surface_roughness_z0;
                // ustar from neighbor cell above
                const amrex::Real ustar =
                    m * kappa /
                    (std::log(1.5_rt * dz / local_z0) - non_neutral_neighbour);
                const amrex::Real ux0 = vel_arrs[nbx](i, j, k, 0);
                const amrex::Real uy0 = vel_arrs[nbx](i, j, k, 1);
                const amrex::Real m0 = std::sqrt((ux0 * ux0) + (uy0 * uy0));
                const amrex::Real uxm1 = vel_arrs[nbx](i, j, k - 1, 0);
                const amrex::Real uym1 = vel_arrs[nbx](i, j, k - 1, 1);
                const amrex::Real mm1 =
                    std::sqrt((uxm1 * uxm1) + (uym1 * uym1));
                const amrex::Real dMdz =
                    amrex::max<amrex::Real>((m0 - mm1) / dz, dMdz_min);
                const amrex::Real mut_loglaw =
                    2.0_rt * ustar * ustar * rho / dMdz;
                const amrex::Real drag =
                    (has_terrain) ? drag_arrs[nbx](i, j, k, 0) : 0.0_rt;
                mu_arrs[nbx](i, j, k) =
                    (mu_arrs[nbx](i, j, k) * (1.0_rt - drag)) +
                    (drag * mut_loglaw);
                const amrex::Real stressScale =
                    (locSurfaceFactor *
                     (std::pow(1.0_rt - fmu, locSurfaceRANSExp) * smag_factor *
                          0.25_rt * locC1 +
                      std::pow(fmu, locSurfaceRANSExp) * ransL)) +
                    ((1.0_rt - locSurfaceFactor) * smag_factor * 0.25_rt *
                     locC1);
                divNij_arrs[nbx](i, j, k, 0) *= rho * stressScale * turnOff *
                                                blankTerrain * non_linear_coeff;
                divNij_arrs[nbx](i, j, k, 1) *= rho * stressScale * turnOff *
                                                blankTerrain * non_linear_coeff;
                divNij_arrs[nbx](i, j, k, 2) *= rho * stressScale * turnOff *
                                                blankTerrain * non_linear_coeff;
            });
    }
    amrex::Gpu::streamSynchronize();

    mu_turb.fillpatch(this->m_sim.time().current_time());
}

/** Kosovic viscosity with the ImmersedTerrain fields.
 *
 *  Same core model as the TerrainDrag kernel, with the terrain fraction and
 *  the surface cells in place of the binary blanking and the drag cells:
 *
 *  - the SGS viscosity and the non-linear term are multiplied by
 *    (1 - w_solid), with w_solid the drag weight of the cell
 *    (``ImmersedTerrain.drag_weight``);
 *  - in surface cells (mask == 2) the viscosity is replaced by the log-law
 *    value 2 rho u*^2 / |dU_t/dn| averaged over the wall patches of the cell
 *    (six faces or surface normal, as ``ImmersedDragForcing.wall_model``),
 *    weighted by (1 - w_solid). The face viscosity is the arithmetic mean of
 *    the cell and the blanked wall-side cell, so this makes the SGS flux
 *    across the interface equal to the wall stress rho u*^2. The friction
 *    velocity uses the same reference cell, distance and normal as the wall
 *    model, and the tangential speed gradient is taken between the cell and
 *    the wall-side cell mirrored from the reference cell.
 */
template <typename Transport>
void Kosovic<Transport>::immersed_terrain_viscosity(
    const int lev,
    const FieldState fstate,
    const ScratchField& gradT,
    const ScratchField& ref_theta)
{
    using immersed_wall::WallPatch;
    using immersedterrain::ImmersedTerrain;

    auto& mu_turb = this->mu_turb();
    const auto& repo = mu_turb.repo();
    const auto& vel = m_vel.state(fstate);
    const auto& den = m_rho.state(fstate);
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> gravity{
        m_gravity[0], m_gravity[1], m_gravity[2]};

    const amrex::Real Cs_sqr = this->m_Cs * this->m_Cs;
    const amrex::Real monin_obukhov_length = m_monin_obukhov_length;
    const amrex::Real dMdz_min = 0.01_rt;
    const amrex::Real locLESTurnOff = m_LESTurnOff;
    const amrex::Real locSwitchLoc = m_switchLoc;
    const amrex::Real locSurfaceRANSExp = m_surfaceRANSExp;
    const amrex::Real locSurfaceFactor = m_surfaceFactor;
    const amrex::Real locC1 = m_C1;
    const amrex::Real tol = constants::TIGHT_TOL;

    const immersed_wall::WallModel wall_model =
        immersed_wall::parse_wall_model(m_ib_wall_model);
    const bool actual_reference = (m_ib_reference_distance == "actual");
    const bool center_weight = (m_ib_drag_weight == "center");
    const amrex::Real solid_threshold = m_ib_solid_threshold;
    const amrex::Real min_z0 = m_ib_min_z0;
    immersed_wall::WallParams wp{};
    wp.kappa = m_kappa;
    wp.beta_m = m_beta_m;
    wp.gamma_m = m_gamma_m;
    wp.inv_mo_length =
        (m_wall_het_model == "mol") ? 1.0_rt / monin_obukhov_length : 0.0_rt;

    const auto& geom = repo.mesh().Geom(lev);
    const auto& problo = geom.ProbLoArray();
    const auto dxv = geom.CellSizeArray();
    const amrex::Real dz = dxv[2];
    const amrex::Real ds = std::cbrt(dxv[0] * dxv[1] * dxv[2]);
    const amrex::Real smag_factor = Cs_sqr * ds * ds;

    const auto& mu_arrs = mu_turb(lev).arrays();
    const auto& rho_arrs = den(lev).const_arrays();
    const auto& vel_arrs = vel(lev).const_arrays();
    const auto& gradT_arrs = gradT(lev).const_arrays();
    const auto& ref_theta_arrs = ref_theta(lev).const_arrays();
    const auto& divNij_arrs = (this->m_divNij)(lev).arrays();
    const auto& frac_arrs =
        this->m_sim.repo().get_field("terrain_fraction")(lev).const_arrays();
    const auto& mask_arrs =
        this->m_sim.repo().get_int_field("terrain_mask")(lev).const_arrays();
    const auto& surf_arrs =
        this->m_sim.repo().get_field("terrain_surface")(lev).const_arrays();
    const auto& z0_arrs =
        this->m_sim.repo().get_field("terrain_roughness")(lev).const_arrays();

    amrex::ParallelFor(
        mu_turb(lev), [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
            const auto& vel_arr = vel_arrs[nbx];
            const auto& frac = frac_arrs[nbx];
            const auto& surf = surf_arrs[nbx];
            const amrex::Real rho = rho_arrs[nbx](i, j, k);
            const amrex::Real beta = frac(i, j, k, 0);
            const amrex::Real w_solid = immersed_wall::solid_weight(
                beta, center_weight, solid_threshold);
            const amrex::Real blankTerrain = 1.0_rt - w_solid;
            const amrex::Real h = surf(i, j, k, ImmersedTerrain::surf_height);
            const amrex::Real z_c = problo[2] + ((k + 0.5_rt) * dz);

            // Kosovic core, as in the TerrainDrag kernel
            const amrex::Real x3 =
                amrex::max<amrex::Real>(z_c - h, 0.5_rt * dz);
            const amrex::Real fmu = std::exp(-x3 / locSwitchLoc);
            const amrex::Real phiM =
                (monin_obukhov_length < 0)
                    ? std::pow(
                          1.0_rt - (16.0_rt * x3 / monin_obukhov_length),
                          -0.25_rt)
                    : 1.0_rt + (5.0_rt * x3 / monin_obukhov_length);
            const amrex::Real wall_distance =
                amrex::max<amrex::Real>(((k + 1) * dz) - h, dz);
            const amrex::Real ransL =
                utils::powi(0.41_rt * wall_distance / phiM, 2);
            const amrex::Real turnOff = std::exp(-x3 / locLESTurnOff);
            const amrex::Real viscosityScale =
                (locSurfaceFactor *
                 (std::pow(1.0_rt - fmu, locSurfaceRANSExp) * smag_factor +
                  std::pow(fmu, locSurfaceRANSExp) * ransL)) +
                ((1.0_rt - locSurfaceFactor) * smag_factor);
            const amrex::Real mut =
                mu_arrs[nbx](i, j, k) * mu_arrs[nbx](i, j, k);
            const amrex::Real T0 = ref_theta_arrs[nbx](i, j, k);
            const amrex::Real stratification_sensor =
                -((gradT_arrs[nbx](i, j, k, 0) * gravity[0]) +
                  (gradT_arrs[nbx](i, j, k, 1) * gravity[1]) +
                  (gradT_arrs[nbx](i, j, k, 2) * gravity[2])) /
                T0;
            const amrex::Real stable_mut =
                mut - (3.0_rt * stratification_sensor);
            amrex::Real stratification = 1.0_rt;
            amrex::Real non_linear_coeff = 1.0_rt;
            if (stratification_sensor > tol) {
                // stable
                non_linear_coeff = (stable_mut < tol) ? 0.0_rt : 1.0_rt;
                stratification =
                    std::sqrt(amrex::max<amrex::Real>(tol, stable_mut));
            } else {
                stratification = std::sqrt(mut);
            }
            mu_arrs[nbx](i, j, k) =
                rho * viscosityScale * turnOff * blankTerrain * stratification;

            // Log-law viscosity on the wall patches of a surface cell
            const bool wall_cell =
                (mask_arrs[nbx](i, j, k, 0) == ImmersedTerrain::mask_surface) &&
                (w_solid < 1.0_rt);
            const amrex::Real z0 =
                amrex::max<amrex::Real>(z0_arrs[nbx](i, j, k, 0), min_z0);
            amrex::GpuArray<WallPatch, 2 * AMREX_SPACEDIM> patches{};
            const int np = wall_cell ? immersed_wall::wall_patches(
                                           wall_model, i, j, k, beta, frac,
                                           surf, dxv, z_c, z0, solid_threshold,
                                           actual_reference, patches.data())
                                     : 0;
            const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> u_c{
                vel_arr(i, j, k, 0), vel_arr(i, j, k, 1), vel_arr(i, j, k, 2)};
            amrex::Real mu_sum = 0.0_rt;
            amrex::Real weight_sum = 0.0_rt;
            for (int ip = 0; ip < np; ++ip) {
                const WallPatch& p = patches[ip];
                const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> u_ref{
                    vel_arr(p.ir, p.jr, p.kr, 0), vel_arr(p.ir, p.jr, p.kr, 1),
                    vel_arr(p.ir, p.jr, p.kr, 2)};
                const amrex::Real ustar = immersed_wall::friction_velocity(
                    immersed_wall::magnitude(
                        immersed_wall::tangential(u_ref, p.nrm)),
                    p, z0, wp);
                // Wall-side cell: mirror of the reference cell
                const int iw = (2 * i) - p.ir;
                const int jw = (2 * j) - p.jr;
                const int kw = (2 * k) - p.kr;
                const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> u_w{
                    vel_arr(iw, jw, kw, 0), vel_arr(iw, jw, kw, 1),
                    vel_arr(iw, jw, kw, 2)};
                const amrex::Real m0 = immersed_wall::magnitude(
                    immersed_wall::tangential(u_c, p.nrm));
                const amrex::Real mw = immersed_wall::magnitude(
                    immersed_wall::tangential(u_w, p.nrm));
                const amrex::Real dMdn =
                    amrex::max<amrex::Real>((m0 - mw) / p.dxn, dMdz_min);
                mu_sum += p.weight * 2.0_rt * ustar * ustar * rho / dMdn;
                weight_sum += p.weight;
            }
            if (weight_sum > 0.0_rt) {
                mu_arrs[nbx](i, j, k) =
                    (1.0_rt - w_solid) * mu_sum / weight_sum;
            }

            const amrex::Real stressScale =
                (locSurfaceFactor *
                 (std::pow(1.0_rt - fmu, locSurfaceRANSExp) * smag_factor *
                      0.25_rt * locC1 +
                  std::pow(fmu, locSurfaceRANSExp) * ransL)) +
                ((1.0_rt - locSurfaceFactor) * smag_factor * 0.25_rt * locC1);
            for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                divNij_arrs[nbx](i, j, k, n) *= rho * stressScale * turnOff *
                                                blankTerrain * non_linear_coeff;
            }
        });
}

template <typename Transport>
void Kosovic<Transport>::update_alphaeff(Field& alphaeff)
{
    BL_PROFILE("kynema-sgf::" + this->identifier() + "::update_alphaeff");

    auto lam_alpha = (this->m_transport).alpha();
    auto& mu_turb = this->m_mu_turb;
    auto& repo = mu_turb.repo();
    const amrex::Real muCoeff = m_muCoeff;
    const int nlevels = repo.num_active_levels();
    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& muturb_arrs = mu_turb(lev).const_arrays();
        const auto& alphaeff_arrs = alphaeff(lev).arrays();
        const auto& lam_diff_arrs = (*lam_alpha)(lev).const_arrays();
        amrex::ParallelFor(
            mu_turb(lev), [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
                alphaeff_arrs[nbx](i, j, k) =
                    lam_diff_arrs[nbx](i, j, k) +
                    (muCoeff * muturb_arrs[nbx](i, j, k));
            });
    }
    amrex::Gpu::streamSynchronize();

    alphaeff.fillpatch(this->m_sim.time().current_time());
}
template <typename Transport>
void Kosovic<Transport>::parse_model_coeffs()
{
    const std::string coeffs_dict = this->model_name() + "_coeffs";
    amrex::ParmParse pp(coeffs_dict);
    pp.query("Cs", this->m_Cs);
}

template <typename Transport>
TurbulenceModel::CoeffsDictType Kosovic<Transport>::model_coeffs() const
{
    return TurbulenceModel::CoeffsDictType{{"Cb", this->m_Cb}};
}

} // namespace turbulence

INSTANTIATE_TURBULENCE_MODEL(Kosovic);

} // namespace kynema_sgf
