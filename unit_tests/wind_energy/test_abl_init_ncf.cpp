#include "abl_test_utils.H"
#include "aw_test_utils/iter_tools.H"
#include "aw_test_utils/test_utils.H"
#include "amr-wind/wind_energy/ABLFieldInit.H"
#include "amr-wind/utilities/ncutils/nc_interface.H"
//#include <filesystem>

namespace amr_wind_tests {
namespace {} // namespace

TEST_F(ABLMeshTest, abl_init_netcdf)
{
    populate_parameters();
    utils::populate_abl_params();

    // Supply name of input NetCDF file to parser
    {
        amrex::ParmParse pp("ABL");
        pp.add("initial_condition_input_file", (std::string) "abl.nc");
    }

    // Create NetCDF file for test
    ncutils::NCFile ncf =
        ncutils::NCFile::create("abl.nc"); //, NC_DISKLESS | NC_NETCDF4);
    // Set up dimensions to correspond to ABL test
    ncf.def_dim("nx", 8);
    ncf.def_dim("ny", 8);
    ncf.def_dim("nz", 64);
    // Define vectors that ABLInit will read
    const std::vector<std::string> three_dim{"nx", "ny", "nz"};
    auto uvel = ncf.def_var("uvel", NC_DOUBLE, three_dim);
    auto vvel = ncf.def_var("vvel", NC_DOUBLE, three_dim);
    auto wvel = ncf.def_var("wvel", NC_DOUBLE, three_dim);
    // Populate std vectors
    const std::vector<size_t> start{0, 0, 0};
    const std::vector<size_t> count{8, 8, 64};
    const std::vector<double> fill_u(8 * 8 * 64, 20.0);
    const std::vector<double> fill_v(8 * 8 * 64, 10.0);
    const std::vector<double> fill_w(8 * 8 * 64, 0.0);
    // Populate NetCDF vectors
    uvel.put(fill_u.data(), start, count);
    vvel.put(fill_v.data(), start, count);
    wvel.put(fill_w.data(), start, count);

    // Set up mock simulation init, starting with mesh and fields
    initialize_mesh();
    auto& frepo = mesh().field_repo();
    auto& velocityf = frepo.declare_field("velocity", 3, 0);
    auto& densityf = frepo.declare_field("density");
    auto& temperaturef = frepo.declare_field("temperature");

    auto velocity = velocityf.vec_ptrs();
    auto density = densityf.vec_ptrs();
    auto temperature = temperaturef.vec_ptrs();

    amr_wind::ABLFieldInit ablinit;
    run_algorithm(
        mesh().num_levels(), density,
        [&](const int lev, const amrex::MFIter& mfi) {
            auto vel = velocity[lev]->array(mfi);
            auto rho = density[lev]->array(mfi);
            auto theta = temperature[lev]->array(mfi);

            const auto& bx = mfi.validbox();
            ablinit(bx, mesh().Geom(lev), vel, rho, theta);
        });

    const int nlevels = mesh().num_levels();
    const amrex::Real tol = 1.0e-12;

    // Test velocity
    {
        amrex::Vector<amrex::Real> min_vel(3), max_vel(3);
        utils::field_minmax(nlevels, velocity, min_vel, max_vel);
        EXPECT_NEAR(min_vel[0], 20.0, tol);
        EXPECT_NEAR(min_vel[1], 10.0, tol);
        EXPECT_NEAR(min_vel[2], 0.0, tol);
        EXPECT_NEAR(max_vel[0], 20.0, tol);
        EXPECT_NEAR(max_vel[1], 10.0, tol);
        EXPECT_NEAR(max_vel[2], 0.0, tol);
    }
}

// Clean up ABL NetCDF file if needed
TEST_F(ABLMeshTest, abl_netcdf_cleanup)
{
    if (std::__fs::filesystem::exists("abl.nc")) {
        remove("abl.nc");
    }
    EXPECT_FALSE(std::__fs::filesystem::exists("abl.nc"));
}

} // namespace amr_wind_tests
