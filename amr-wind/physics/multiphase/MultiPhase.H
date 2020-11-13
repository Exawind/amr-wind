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

enum class multiphase_model { VOF, LS, CLSVOF };

class MultiPhase : public Physics::Register<MultiPhase>
{
public:
    static const std::string identifier() { return "MultiPhase"; }

    MultiPhase(CFDSim& sim);

    virtual ~MultiPhase() = default;

    void initialize_fields(int, const amrex::Geometry&) override {}

    void post_init_actions() override;

    void post_regrid_actions() override {}

    void pre_advance_work() override {}

    void post_advance_work() override;

    void set_multiphase_properties(int level, const amrex::Geometry& geom);

private:
    const CFDSim& m_sim;

    Field& m_velocity;
    Field& m_mueff;
    Field& m_density;

    // Pointer to the levelset
    Field* m_levelset{nullptr};

    // Pointer to VOF
    Field* m_vof{nullptr};

    // Density value for Fluid 1
    amrex::Real m_rho1{10.0};

    // Density value for Fluid 2
    amrex::Real m_rho2{1.0};

    // Dynamic Viscosity for Fluid 1
    amrex::Real m_mu1{1.0};

    // Dynamic Viscosity for Fluid 2
    amrex::Real m_mu2{1.0};

    // Multiphase model enum
    multiphase_model m_multiphase_model=multiphase_model::VOF;
};

} // namespace amr_wind

#endif /* MULTIPHASE_H */
