#include "gtest/gtest.h"
#include "aw_test_utils/MeshTest.H"

#include "PDETraits.H"
#include "SchemeTraits.H"
#include "PDE.H"
#include "icns/icns_ops.H"

namespace amr_wind_tests {

class PDETest : public MeshTest
{};

TEST_F(PDETest, test_pde_impl)
{
    initialize_mesh();

    amr_wind::SimTime time;

    using ICNSGodunov =
        amr_wind::pde::PDESystem<amr_wind::pde::ICNS, amr_wind::fvm::Godunov>;
    using TemperatureGodunov =
        amr_wind::pde::PDESystem<amr_wind::pde::Temperature, amr_wind::fvm::Godunov>;

    ICNSGodunov lowmach(time, mesh().field_repo(), 35);
    TemperatureGodunov theta(time, mesh().field_repo(), 35);

    amrex::Print()
        << "Num registered fields: "
        << mesh().field_repo().num_fields() << std::endl;

    for (auto& ff: mesh().field_repo().fields()) {
        amrex::Print() << " - " << ff->name() << std::endl;
    }

    lowmach.compute_source_term(amr_wind::FieldState::Old);
    theta.compute_source_term(amr_wind::FieldState::Old);
}

}
