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

#include "amr-wind/utilities/FieldPlaneAveraging.H"
#include "amr-wind/utilities/trig_ops.H"

namespace amr_wind {

namespace {

auto reflect(vs::Vector line, vs::Vector vec)
{

    vs::Tensor ref(
        1 - 2 * line.x() * line.x(), -2 * line.x() * line.y(),
        -2 * line.x() * line.z(), -2 * line.y() * line.x(),
        1 - 2 * line.y() * line.y(), -2 * line.y() * line.z(),
        -2 * line.z() * line.x(), -2 * line.z() * line.y(),
        1 - 2 * line.z() * line.z());

    return vec & ref;
}

auto rotate_euler_vec(vs::Vector axis, double angle, vs::Vector vec)
{

    axis.normalize();

    const auto RotMat = vs::quaternion(axis, angle);

    return vec & RotMat;
}

vs::Vector rotation(const vs::Vector& angles, const vs::Vector& data)
{
    const vs::Tensor rotMatrix =
        vs::xrot(angles.x()) & vs::yrot(angles.y()) & vs::zrot(angles.z());
    return data & rotMatrix;
}

auto AdjustLidarPattern(
    vs::Vector beamPt,
    double yaw,
    double pitch,
    double roll,
    vs::Vector translation)
{

    const vs::Vector angles(pitch, yaw, roll);

    auto beamPt_transform = rotation(angles, beamPt) + translation;

    return beamPt_transform;
}

auto generate_pattern(double time)
{
    vs::Vector axis(1, 0, 0);

    const double innerPrism_theta0 = 90 * M_PI / 180;
    const double outerPrism_theta0 = 90 * M_PI / 180;
    const double innerPrism_rot = 3.5;
    const double outerPrism_rot = 6.5;
    const double innerPrism_azimuth = 15.2 * M_PI / 180;
    const double outerPrism_azimuth = 15.2 * M_PI / 180;
    const vs::Vector ground(0, 0, 1);

    const double innerTheta =
        innerPrism_theta0 + innerPrism_rot * time * 2 * M_PI;
    const double outerTheta =
        outerPrism_theta0 + outerPrism_rot * time * 2 * M_PI;

    const auto reflection_1 = rotate_euler_vec(
        axis, innerTheta,
        rotate_euler_vec(ground, -(innerPrism_azimuth / 2 + M_PI / 2), axis));

    const auto reflection_2 = rotate_euler_vec(
        axis, outerTheta,
        rotate_euler_vec(ground, outerPrism_azimuth / 2, axis));

    return reflect(reflection_2, reflect(reflection_1, axis));
}
} // namespace
TEST(SpinnerTest, test_spinner)
{
    double time = 0;
    double yaw = 180;
    double pitch = 0;
    double roll = 0;
    vs::Vector translation(0, 0, 0);
    std::ofstream outputFile("SpinnerLidarpattern.txt");
    outputFile << "x,y,z" << std::endl;

    for (int j = 0; j < 100; ++j) {
        auto beam_vector = generate_pattern(2. / 100 * j);
        beam_vector =
            AdjustLidarPattern(beam_vector, pitch, roll, yaw, translation);
        outputFile << beam_vector[0] << "," << beam_vector[1] << ","
                   << beam_vector[2] << std::endl;
    }
}
} // namespace amr_wind
