#include "aw_test_utils/MeshTest.H"
#include "aw_test_utils/iter_tools.H"

#include "AMReX_Box.H"
#include "AMReX_BoxArray.H"
#include "AMReX_BoxList.H"
#include "AMReX_Geometry.H"
#include "AMReX_RealBox.H"
#include "AMReX_Vector.H"
#include "amr-wind/core/vs/tensorI.H"
#include "amr-wind/core/vs/tensor.H"
#include "amr-wind/utilities/sampling/DTUSpinnerSampler.H"

#include "amr-wind/utilities/FieldPlaneAveraging.H"
#include "amr-wind/utilities/trig_ops.H"

namespace amr_wind::sampling {

TEST(SpinnerTest, test_spinner)
{
    double yaw = 180;
    double pitch = 0;
    double roll = 0;
    vs::Vector translation(0, 0, 0);
    std::ofstream outputFile("SpinnerLidarpattern.txt");
    outputFile << "x,y,z" << std::endl;

    PrismParameters InnerPrism;
    PrismParameters OuterPrism;
    OuterPrism.rot = 6.5;

    for (int j = 0; j < 983; ++j) {
        auto beam_vector =
            generate_lidar_pattern(InnerPrism, OuterPrism, 2. / 984 * j);
        beam_vector = adjust_lidar_pattern(beam_vector, yaw, pitch, roll);
        outputFile << beam_vector[0] << "," << beam_vector[1] << ","
                   << beam_vector[2] << std::endl;
    }
}

} // namespace amr_wind::sampling
