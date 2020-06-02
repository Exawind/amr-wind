#include "physics_test_utils.H"

namespace amr_wind_tests {

PhysicsEx::PhysicsEx(amr_wind::CFDSim&) {}

void PhysicsEx::post_init_actions()  {}
void PhysicsEx::post_regrid_actions()  {}
void PhysicsEx::initialize_fields(int, const amrex::Geometry&)  {}
void PhysicsEx::pre_advance_work()  {}
void PhysicsEx::post_advance_work()  {}

}
