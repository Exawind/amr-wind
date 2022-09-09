#include "ocean_waves_test_utils.H"
#include "amr-wind/utilities/trig_ops.H"
#include "aw_test_utils/iter_tools.H"
#include "aw_test_utils/test_utils.H"

namespace amr_wind_tests {

TEST_F(OceanWavesMeshTest, WaveGeneration)
{
    constexpr amrex::Real tol = 1.0e-12;
    populate_parameters();
    utils::populate_ocean_waves_params();

    initialize_mesh();
    auto& frepo = mesh().field_repo();
    auto& velocityf = frepo.declare_field("velocity", 3, 0);
    auto& densityf = frepo.declare_field("density");
    auto& voff = frepo.declare_field("vof");
}

TEST_F(OceanWavesMeshTest, WaveAbsorption)
{
    constexpr amrex::Real tol = 1.0e-12;
    populate_parameters();
    utils::populate_ocean_waves_params();

    initialize_mesh();
    auto& frepo = mesh().field_repo();
    auto& velocityf = frepo.declare_field("velocity", 3, 0);
    auto& densityf = frepo.declare_field("density");
    auto& voff = frepo.declare_field("vof");
}

} // namespace amr_wind_tests
