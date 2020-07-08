#ifndef ABL_TEST_UTILS_H
#define ABL_TEST_UTILS_H

#include "aw_test_utils/MeshTest.H"
#include "aw_test_utils/pp_utils.H"

namespace amr_wind_tests {

class ABLTest : public AmrexTest
{};

// Introduce namespace for the text fixture
class ABLMeshTest : public MeshTest
{};

namespace utils {
void populate_abl_params();
}
} // namespace amr_wind_tests

#endif /* ABL_TEST_UTILS_H */
