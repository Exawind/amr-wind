#ifndef KLAXELL_H
#define KLAXELL_H

#include <string>
#include "src/core/ScratchField.H"
#include "src/turbulence/TurbModelBase.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace kynema_sgf::turbulence {

/** Base class for 1-Equation RANS TKE turbulence model
 *  \ingroup turb_model
 *
 *  With terrain the mixing length uses the height above the terrain and the
 *  eddy viscosity is zero inside the terrain. ``KLAxell.terrain_model``
 *  selects the terrain physics this is evaluated from (default
 *  ``TerrainDrag``; ``ImmersedTerrain`` weights the viscosity by one minus
 *  the drag weight of the cell in a separate kernel).
 */
template <typename Transport>
class KLAxell : public TurbModelBase<Transport>
{
public:
    static std::string identifier()
    {
        return "KLAxell-" + Transport::identifier();
    }

    explicit KLAxell(CFDSim& sim);

    [[nodiscard]] std::string model_name() const override { return "KLAxell"; }

    //! Update the turbulent viscosity field
    void update_turbulent_viscosity(
        FieldState fstate, DiffusionType /*unused*/) override;

    //! Do any post advance work
    void post_advance_work() override;

    //! Actions to perform post initialization
    void post_init_actions() override;

    //! Actions to perform post regrid
    void post_regrid_actions() override;

    //! Update the effective thermal diffusivity field
    void update_alphaeff(Field& alphaeff) override;

    //! Turbulent viscosity on one level from the ImmersedTerrain fields,
    //! used when terrain_model is ImmersedTerrain
    void immersed_terrain_viscosity(
        int lev,
        FieldState fstate,
        const ScratchField& gradT,
        const ScratchField& beta_field);

    //! Update the effective scalar diffusivity field
    void update_scalar_diff(Field& deff, const std::string& name) override;

    //! Parse turbulence model coefficients
    void parse_model_coeffs() override;

    //! Return turbulence model coefficients
    [[nodiscard]] TurbulenceModel::CoeffsDictType model_coeffs() const override;

private:
    Field& m_vel;
    Field& m_turb_lscale;
    Field& m_shear_prod;
    Field& m_buoy_prod;
    Field& m_dissip;
    Field& m_rho;
    Field* m_tke{nullptr};
    //! Turbulence constant
    amrex::Real m_Cmu{0.556_rt};
    amrex::Real m_Cmu_prime{0.556_rt};
    amrex::Real m_Cb_stable{0.25_rt};
    amrex::Real m_Cb_unstable{0.35_rt};
    amrex::Real m_prandtl{1.0_rt};
    Field& m_temperature;
    //! Gravity vector (m/s^2)
    amrex::Vector<amrex::Real> m_gravity{0.0_rt, 0.0_rt, -9.81_rt};
    amrex::Real m_surf_flux{0};
    amrex::Real m_meso_sponge_start{2000};

    std::unique_ptr<ScratchField> m_gradT;

    //! Terrain physics the model couples to: ``TerrainDrag`` (default, binary
    //! blanking) or ``ImmersedTerrain`` (partial fraction). Also read by the
    //! ``KransAxell`` TKE source.
    std::string m_terrain_model{"TerrainDrag"};
    //! Surface-cell settings shared with ImmersedTerrain
    amrex::Real m_ib_solid_threshold{0.5_rt};
    std::string m_ib_drag_weight{"fraction"};
};

} // namespace kynema_sgf::turbulence

#endif /* KLAXELL_H */
