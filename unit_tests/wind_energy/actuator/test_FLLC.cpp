#include "aw_test_utils/AmrexTest.H"
#include "amr-wind/wind_energy/actuator/FLLC.H"
#include "amr-wind/wind_energy/actuator/FLLCOpBase.H"

namespace amr_wind {
namespace actuator {

TEST(TestFLLCData, data_sized_by_single_function)
{
    const int num_points = 10;

    VecList dummy_vec(num_points);
    RealList dummy_real(num_points);
    TensorList dummy_tensor(num_points);

    FLLCData data(dummy_vec);

    ComponentView view;
    view.pos = ::amr_wind::utils::slice(dummy_vec, 0, num_points);
    view.force = ::amr_wind::utils::slice(dummy_vec, 0, num_points);
    view.epsilon = ::amr_wind::utils::slice(dummy_vec, 0, num_points);
    view.orientation = ::amr_wind::utils::slice(dummy_tensor, 0, num_points);
    view.chord = ::amr_wind::utils::slice(dummy_real, 0, num_points);

    data.init_data(view, 1.0);

    ASSERT_EQ(num_points, data.les_velocity.size());
    ASSERT_EQ(num_points, data.optimal_velocity.size());
    ASSERT_EQ(num_points, data.lift.size());
    ASSERT_EQ(num_points, data.grad_lift.size());
}

} // namespace actuator
} // namespace amr_wind
