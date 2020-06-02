#include "gtest/gtest.h"
#include "aw_test_utils/MeshTest.H"
#include "amr-wind/turbulence/TurbulenceModel.H"

namespace amr_wind_tests {

class TurbTest : public MeshTest
{};

TEST_F(TurbTest, test_turb_create)
{
    initialize_mesh();

    amr_wind::turbulence::TurbulenceModel::print(std::cout);
}

}
