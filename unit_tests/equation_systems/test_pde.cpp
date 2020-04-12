#include "gtest/gtest.h"
#include "aw_test_utils/MeshTest.H"

#include "PDETraits.H"
#include "SchemeTraits.H"
#include "PDE.H"
#include "icns/icns_ops.H"

namespace amr_wind_tests {

class PDETest : public MeshTest
{};

TEST_F(PDETest, test_pde_create_godunov)
{
    initialize_mesh();
    amr_wind::SimTime time;

    auto lowmach = amr_wind::pde::PDEBase::create(
        "ICNS-Godunov", time, mesh().field_repo(), 35);
    auto theta = amr_wind::pde::PDEBase::create(
        "Temperature-Godunov", time, mesh().field_repo(), 35);

    EXPECT_EQ(mesh().field_repo().num_fields(), 19);
}

TEST_F(PDETest, test_pde_create_mol)
{
    initialize_mesh();
    amr_wind::SimTime time;

    auto lowmach = amr_wind::pde::PDEBase::create(
        "ICNS-MOL", time, mesh().field_repo(), 35);
    auto theta = amr_wind::pde::PDEBase::create(
        "Temperature-MOL", time, mesh().field_repo(), 35);

    EXPECT_EQ(mesh().field_repo().num_fields(), 23);
}

}
