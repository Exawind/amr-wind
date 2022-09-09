#include "ocean_waves_test_utils.H"
#include "amr-wind/utilities/trig_ops.H"
#include "aw_test_utils/iter_tools.H"
#include "aw_test_utils/test_utils.H"

namespace amr_wind_tests {

TEST_F(OceanWavesMeshTest, WaveGeneration)
{
    constexpr amrex::Real tol = 1.0e-12;
    utils::populate_ocean_waves_params();
}

} // namespace amr_wind_tests
