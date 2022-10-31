#ifndef MULTIPHASE_H
#define MULTIPHASE_H

#include "amr-wind/core/Physics.H"
#include "amr-wind/core/Field.H"

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

    void post_regrid_actions() override {}

    void pre_advance_work() override;

    void post_advance_work() override;

    void set_density_via_levelset();

    void set_density_via_vof();

    void set_nph_density();

    void calculate_advected_facedensity();

    void levelset2vof();

    void favre_filtering();

    amrex::Real volume_fraction_sum();

    amrex::Real momentum_sum(int n);

    InterfaceCapturingMethod interface_capturing_method();

    amrex::Real rho1() const { return m_rho1; }

    amrex::Real rho2() const { return m_rho2; }

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

    bool m_interface_smoothing{false};

    int m_smooth_freq{10};

    // Multiphase model enum
    std::string m_interface_model{"vof"};
    InterfaceCapturingMethod m_interface_capturing_method =
        InterfaceCapturingMethod::VOF;

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
