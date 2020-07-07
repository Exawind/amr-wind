#include "gtest/gtest.h"
#include "amr-wind/utilities/ncutils/nc_interface.H"
#include "AMReX.H"

namespace amr_wind_tests {

TEST(NetCDFUtils, ncfile)
{
    constexpr int num_points = 10;
    ncutils::NCFile ncf =
        ncutils::NCFile::create("test.nc", NC_DISKLESS | NC_NETCDF4);
    ASSERT_GE(ncf.ncid, 0);

    auto tdim = ncf.def_dim("nsteps", NC_UNLIMITED);
    auto dim1 = ncf.def_dim("nx", 10);
    auto dim2 = ncf.dim("nx");
    auto var1 = ncf.def_array("vel", NC_DOUBLE, {"nx"});
    auto var2 = ncf.def_array("timehist", NC_DOUBLE, {"nsteps", "nx"});

    {
        ASSERT_EQ(dim1.dimid, dim2.dimid);
        ASSERT_EQ(dim1.len(), num_points);
        ASSERT_EQ(var1.ndim(), 1);

        ASSERT_EQ(ncf.num_dimensions(), 2);
        ASSERT_EQ(ncf.num_variables(), 2);

        auto vshp = var1.shape();
        ASSERT_EQ(vshp.size(), 1u);
        ASSERT_EQ(vshp[0], num_points);

        auto vshp2 = var2.shape();
        ASSERT_EQ(vshp2.size(), 2u);
        ASSERT_EQ(vshp2[0], 0);
        ASSERT_EQ(vshp2[1], num_points);
    }

    std::vector<double> fill_val(num_points, 2000.0);
    std::vector<size_t> start{0, 0};
    std::vector<size_t> count{1, num_points};

    var1.put(fill_val.data());
    var2.put(fill_val.data(), start, count);

    {
        ASSERT_EQ(tdim.len(), 1u);
        auto vshp2 = var2.shape();
        ASSERT_EQ(vshp2.size(), 2u);
        ASSERT_EQ(vshp2[0], 1);
        ASSERT_EQ(vshp2[1], num_points);
    }
}

TEST(NetCDFUtils, ncgroups)
{
    constexpr int num_points = 10;
    ncutils::NCFile ncf =
        ncutils::NCFile::create("test_groups.nc", NC_DISKLESS | NC_NETCDF4);
    ASSERT_GE(ncf.ncid, 0);

    ncf.def_dim("nsteps", NC_UNLIMITED);
    ncf.def_dim("nx", num_points);
    auto line1 = ncf.def_group("line1");
    line1.def_dim("ny", num_points);
    auto vel = line1.def_var("vel", NC_DOUBLE, {"nx", "ny"});
    std::vector<double> fill_val(num_points * num_points, 2000.0);
    vel.put(fill_val.data());

    ASSERT_EQ(ncf.num_groups(), 1);
    ASSERT_EQ(ncf.num_dimensions(), 2);
    ASSERT_EQ(line1.num_dimensions(), 1);

    {
        auto vshape = vel.shape();
        ASSERT_EQ(vshape.size(), 2u);
        ASSERT_EQ(vshape[0], num_points);
        ASSERT_EQ(vshape[1], num_points);
    }

    ASSERT_TRUE(ncf.has_group("line1"));
    ASSERT_FALSE(ncf.has_group("line2"));

    // The file is closed when the NCFile variable goes out of scope. This test
    // checks if multiple calls to close are safe.
    ncf.close();

    // These tests require an actual file on disk... doesn't work with diskless mode
#if 0
    {
        auto ncf1 = ncutils::NCFile::open("test_groups.nc", NC_WRITE);
        ASSERT_EQ(ncf1.num_groups(), 1);
        ASSERT_EQ(ncf1.num_dimensions(), 2);
        auto line1 = ncf1.group("line1");
        ASSERT_EQ(line1.num_variables(), 1);
        auto vel = line1.var("vel");
        ASSERT_EQ(vel.ndim(), 2u);
    }
#endif
}

TEST(NetCDFUtils, var_io)
{
    constexpr int num_points = 5;
    ncutils::NCFile ncf =
        ncutils::NCFile::create("test_vario.nc", NC_DISKLESS | NC_NETCDF4);
    ASSERT_GE(ncf.ncid, 0);
    ncf.put_attr("title", "AMR-Wind data sampling");

    auto line1 = ncf.def_group("line1");
    line1.put_attr("sampling_type", "LineSampler");
    line1.def_dim("nx", num_points);
    line1.def_dim("ny", num_points);
    std::vector<double> fill_val(num_points * num_points);
    ASSERT_EQ(line1.get_attr("sampling_type"), "LineSampler");

    int idx = 0;
    for (int i=0; i < num_points; ++i)
        for (int j=0; j < num_points; ++j)
            fill_val[idx++] = static_cast<double>(i * num_points + j);

    auto vel = line1.def_var("vel", NC_DOUBLE, {"nx", "ny"});
    vel.put_attr("units", "m/s");
    vel.put(fill_val.data());
    ASSERT_EQ("m/s", vel.get_attr("units"));

    const size_t istart = num_points - 2;
    const double start_val = static_cast<double>(istart * num_points);
    std::vector<double> buf(num_points);
    vel.get(buf.data(), {istart, 0}, {1, num_points});

    for (int i=0; i < num_points; ++i)
        ASSERT_NEAR(buf[i], start_val + static_cast<double>(i), 1.0e-12);

    // Test hyperslab
    buf.resize(num_points*2);
    vel.get(buf.data(), {0, 0}, {num_points, 2}, {1, 2});
    idx = 0;
    for (int i=0; i < num_points; ++i)
        for (int j=0; j < 2; ++j)
            ASSERT_NEAR(buf[idx++], static_cast<double>(num_points * i + 2 * j), 1.0e-12);
}

} // namespace amr_wind_tests
