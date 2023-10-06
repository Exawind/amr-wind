#ifndef MULTIPHASE_H
#define MULTIPHASE_H

#include "amr-wind/core/Physics.H"
#include "amr-wind/core/Field.H"
#include "amr-wind/core/IntField.H"
#include "amr-wind/core/ScratchField.H"

/** Multiphase physics
 *
 * This class implemnts the necessary code to simulate multiphase
 * (air-water interaction) problems with AMR-Wind. It reads the
 * namespace `MultiPhase' from the input file to determine how the
 * initial conditions are generated.
 *
 */
namespace amr_wind {

enum class InterfaceCapturingMethod {
    VOF, ///< Volume of fluid
    LS   ///< Levelset
};

class MultiPhase : public Physics::Register<MultiPhase>
{
public:
    static std::string identifier() { return "MultiPhase"; }

    explicit MultiPhase(CFDSim& sim);

    ~MultiPhase() override = default;

    void
    initialize_fields(int /*level*/, const amrex::Geometry& /*geom*/) override
    {}

    void post_init_actions() override;

    void post_regrid_actions() override;

    void pre_advance_work() override;

    void post_advance_work() override;

    void set_density_via_levelset();

    void set_density_via_vof(
        amr_wind::FieldState fstate = amr_wind::FieldState::New);

    void set_nph_density();

    void calculate_advected_facedensity();

    void levelset2vof();

    void levelset2vof(IntField& iblank_cell, ScratchField& vof_scr);

    void favre_filtering();

    void set_intf_init_method_to_vof()
    {
        m_interface_init_method = InterfaceCapturingMethod::VOF;
    }

    void set_intf_init_method_to_ls()
    {
        m_interface_init_method = InterfaceCapturingMethod::LS;
    };

    amrex::Real volume_fraction_sum();

    amrex::Real momentum_sum(int n);

    InterfaceCapturingMethod interface_capturing_method();

    amrex::Real rho1() const { return m_rho1; }

    amrex::Real rho2() const { return m_rho2; }

    amrex::Vector<amrex::Real> gravity() const { return m_gravity; }

private:
    const CFDSim& m_sim;

    Field& m_velocity;
    Field& m_density;

    // Pointer to the levelset
    Field* m_levelset{nullptr};

    // Pointer to VOF
    Field* m_vof{nullptr};

    // Density value for Fluid 1
    amrex::Real m_rho1{10.0};

    // Density value for Fluid 2
    amrex::Real m_rho2{1.0};

    // Bools according to ICNS settings
    // Turning on perturbational density and pressure
    bool is_pptb{false};
    // Reconstructing true pressure field at end of timestep
    bool is_ptrue{false};

    bool m_init_p{false};

    // Info to create rho0
    amrex::Real water_level0{0.0};
    // Info to reconstruct true pressure
    amrex::Vector<amrex::Real> m_gravity{{0.0, 0.0, -9.81}};

    bool m_interface_smoothing{false};

    int m_smooth_freq{10};

    // Multiphase model enum
    std::string m_interface_model{"vof"};
    InterfaceCapturingMethod m_interface_capturing_method =
        InterfaceCapturingMethod::VOF;

    InterfaceCapturingMethod m_interface_init_method =
        InterfaceCapturingMethod::LS;

    // Verbose flag for multiphase
    int m_verbose{0};

    // sum of volume fractions (for vof only)
    amrex::Real m_total_volfrac{0.0};

    amrex::Real q0{0.0};
    amrex::Real q1{0.0};
    amrex::Real q2{0.0};
    amrex::Real sumvof0{0.0};
};

} // namespace amr_wind

#endif /* MULTIPHASE_H */
