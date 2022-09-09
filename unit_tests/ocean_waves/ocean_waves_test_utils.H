#ifndef OCEAN_WAVES_TEST_UTILS_H
#define OCEAN_WAVES_TEST_UTILS_H

#include "aw_test_utils/MeshTest.H"
#include "aw_test_utils/pp_utils.H"

namespace amr_wind_tests {

class OceanWavesTest : public AmrexTest
{};

// Introduce namespace for the text fixture
class OceanWavesMeshTest : public MeshTest
{};

namespace utils {
void populate_ocean_waves_params();
}
} // namespace amr_wind_tests

#endif /* OCEAN_WAVES_TEST_UTILS_H */
