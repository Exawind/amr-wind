#ifndef FuzzyInterface_H
#define FuzzyInterface_H

#include "amr-wind/core/Physics.H"
#include "amr-wind/core/Field.H"

/** Multiphase 3D physics
 *  \ingroup multiphase_physics
 *
 */

namespace amr_wind {

class FuzzyInterface : public Physics::Register<FuzzyInterface>
{
public:
    static std::string identifier() { return "FuzzyInterface"; }

    explicit FuzzyInterface(CFDSim& sim);

    ~FuzzyInterface() override = default;

    //! Initialize the levelset and velocity fields for Sloshing Tank
    //! simulations
    void initialize_fields(int level, const amrex::Geometry& geom) override;

    void post_init_actions() override {}

    void post_regrid_actions() override {}

    void pre_advance_work() override {}

    void post_advance_work() override {}

private:
    Field& m_velocity;
    Field& m_vof;
    Field& m_pressure;

    //! Initial zero-level free-surface water depth
    amrex::Real m_waterlevel{0.0};
    //! Divide vof between smooth and discrete at x
    amrex::Real m_lx_vofjump{-0.25};
    amrex::Real m_hx_vofjump{0.25};
    // Thickness of interface in smooth section
    amrex::Real m_intf_th{0.0};

    //! Stuff to get from MultiPhase physics
    amrex::Real m_rho1{1000.};
    amrex::Real m_rho2{1.};
    amrex::Vector<amrex::Real> m_gravity{0.0, 0.0, -9.81};
};

} // namespace amr_wind
#endif /* FuzzyInterface_H */
