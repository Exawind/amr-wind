#include "aw_test_utils/MeshTest.H"
#include "amr-wind/fvm/gradient.H"

namespace amr_wind_tests {

class FvmOpTest : public MeshTest
{};

TEST_F(FvmOpTest, gradient)
{
    initialize_mesh();

    auto& repo = sim().repo();
    auto& vel = repo.declare_field("vel", 3, 1);

    auto grad_vel = amr_wind::fvm::gradient(vel);
}

}
