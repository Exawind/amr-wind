#ifndef PHYSICS_TEST_UTILS_H
#define PHYSICS_TEST_UTILS_H

#include <string>

#include "amr-wind/core/Physics.H"

namespace amr_wind_tests {

class PhysicsEx : public amr_wind::Physics::Register<PhysicsEx>
{
public:
    static const std::string identifier() { return "PhysicsEx"; }

    PhysicsEx(amr_wind::CFDSim&);

    virtual ~PhysicsEx() = default;

    void post_init_actions() override;
    void post_regrid_actions() override;
    void initialize_fields(int, const amrex::Geometry&) override;
    void pre_advance_work() override;
    void post_advance_work() override;
};

}

#endif /* PHYSICS_TEST_UTILS_H */
