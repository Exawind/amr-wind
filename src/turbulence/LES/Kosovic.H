#ifndef KOSOVIC_H
#define KOSOVIC_H

#include <string>
#include "src/turbulence/TurbModelBase.H"
#include "src/core/FieldRepo.H"
#include "src/utilities/constants.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace kynema_sgf::turbulence {

/** Kosovic LES Model
 *  \ingroup turb_model
 *
 *  With terrain the SGS viscosity is zero inside the terrain and, in the
 *  first fluid cells, replaced by the value that makes the SGS flux across
 *  the terrain interface equal to the log-law wall stress rho u*^2 (the face
 *  viscosity is the arithmetic mean of the cell and the blanked neighbor, so
 *  the cell value is 2 rho u*^2 / |dU/dn|). ``Kosovic.terrain_model``
 *  selects the terrain physics this is evaluated from.
 */
template <typename Transport>
class Kosovic : public TurbModelBase<Transport>
{
public:
    static std::string identifier()
    {
        return "Kosovic-" + Transport::identifier();
    }

    explicit Kosovic(CFDSim& sim);

    //! Model name for debugging purposes
    [[nodiscard]] std::string model_name() const override { return "Kosovic"; }

    //! Update the turbulent viscosity field
    void update_turbulent_viscosity(
        FieldState fstate, DiffusionType /*unused*/) override;

    //! No post advance work for this model
    void post_advance_work() override {}

    //! Update the effective thermal diffusivity field
    void update_alphaeff(Field& alphaeff) override;

    //! Turbulent viscosity on one level with the ImmersedTerrain fields
    //! (``terrain_model = ImmersedTerrain``)
    void immersed_terrain_viscosity(
        int lev,
        FieldState fstate,
        const ScratchField& gradT,
        const ScratchField& ref_theta);

    //! Parse turbulence model coefficients
    void parse_model_coeffs() override;

    //! Return model coefficients dictionary
    [[nodiscard]] TurbulenceModel::CoeffsDictType model_coeffs() const override;

private:
    //! Kosovic coefficient (default value set for ABL simulations)
    // Ref: Mirocha et. al "Implementation of a Nonlinear Subfilter Turbulence
    //  Stress Model for Large-Eddy Simulation in the Advanced Research WRF
    //  Model"
    // , MWR 2012.
    amrex::Real m_Cb{0.36_rt};
    amrex::Real m_Cs{0.135_rt};
    amrex::Real m_C1{2.1_rt};
    amrex::Real m_C2{2.1_rt};
    amrex::Real m_Sk{0.5_rt};
    bool m_surfaceRANS{false};
    amrex::Real m_surfaceFactor{0.0_rt};
    amrex::Real m_switchLoc{24.0_rt};
    amrex::Real m_surfaceRANSExp{2.0_rt};
    amrex::Real m_LESTurnOff{constants::LARGE_NUM};
    amrex::Real m_muCoeff{3.0_rt};
    bool m_writeTerms{false};
    const Field& m_vel;
    const Field& m_rho;
    const Field& m_theta;
    Field& m_Nij;
    Field& m_divNij;
    std::string m_wall_het_model{"none"};
    amrex::Real m_monin_obukhov_length{constants::LARGE_NUM};
    amrex::Real m_kappa{0.41_rt};
    amrex::Real m_gamma_m{5.0_rt};
    amrex::Real m_beta_m{16.0_rt};
    amrex::Real m_surface_roughness_z0{0.1_rt};
    amrex::Vector<amrex::Real> m_gravity{0.0_rt, 0.0_rt, -9.81_rt};

    //! Terrain physics the SGS model couples to: ``TerrainDrag`` (default,
    //! binary blanking, log-law viscosity in the drag cells above the
    //! terrain) or ``ImmersedTerrain`` (partial fraction, log-law viscosity
    //! on the wall patches of the surface cells, six faces or surface
    //! normal, as in ImmersedDragForcing)
    std::string m_terrain_model{"TerrainDrag"};
    //! Wall-model settings shared with ImmersedDragForcing
    std::string m_ib_wall_model{"cell_offset"};
    std::string m_ib_reference_distance{"nominal"};
    amrex::Real m_ib_min_z0{1.0e-4_rt};
    //! Surface-cell settings shared with ImmersedTerrain
    amrex::Real m_ib_solid_threshold{0.5_rt};
    std::string m_ib_drag_weight{"fraction"};
};

} // namespace kynema_sgf::turbulence

#endif /* SMAGORINSKY_H */
