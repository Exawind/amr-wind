
#include "aw_test_utils/MeshTest.H"

#include "amr-wind/utilities/sampling/FreeSurface.H"

namespace amr_wind_tests {

namespace {

amrex::Real init_vof(amr_wind::Field& vof_fld, amrex::Real water_level)
{
    const auto& mesh = vof_fld.repo().mesh();
    const int nlevels = vof_fld.repo().num_active_levels();

    // Since VOF is cell centered
    amrex::Real offset = 0.5;
    // VOF has only one component
    int d = 0;
    // Linearly interpolated water level
    amrex::Real liwl_max = -100.0;

    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& dx = mesh.Geom(lev).CellSizeArray();
        const auto& problo = mesh.Geom(lev).ProbLoArray();

        for (amrex::MFIter mfi(vof_fld(lev)); mfi.isValid(); ++mfi) {
            auto bx = mfi.growntilebox();
            const auto& farr = vof_fld(lev).array(mfi);

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                const amrex::Real z = problo[2] + (k + offset) * dx[2];

                amrex::Real local_vof = std::min(
                    1.0,
                    std::max(
                        0.0, (water_level - (z - offset * dx[2])) / dx[2]));
                farr(i, j, k, d) = local_vof;
            });
        }
        // Use ReduceMax to get linearly interpolated water level
        amrex::Real liwl_lev = -100.0;
        liwl_lev = amrex::ReduceMax(
            vof_fld(lev), 0.,
            [=] AMREX_GPU_HOST_DEVICE(
                amrex::Box const& b,
                amrex::Array4<amrex::Real const> const& field_arr)
                -> amrex::Real {
                amrex::Real liwl = -100.0;
                amrex::Loop(b, [=, &liwl](int i, int j, int k) noexcept {
                    const amrex::Real z = problo[2] + (k + offset) * dx[2];
                    auto local_vof = field_arr(i, j, k, 0);
                    if (local_vof > 0 && local_vof < 1) {
                        liwl =
                            z + (0.5 - local_vof) / (0.0 - local_vof) * dx[2];
                    }
                });
                return liwl;
            });
        liwl_max = amrex::max(liwl_max, liwl_lev);
    }
    amrex::ParallelDescriptor::ReduceRealMax(liwl_max);
    // If only single-phase cells, no interpolation needed
    if (liwl_max < 0) liwl_max = water_level;
    // Return result
    return liwl_max;
}

class FreeSurfaceImpl : public amr_wind::free_surface::FreeSurface
{
public:
    FreeSurfaceImpl(amr_wind::CFDSim& sim, const std::string& label)
        : amr_wind::free_surface::FreeSurface(sim, label)
    {}
    int check_output(const std::string& op, amrex::Real check_val);
    int check_pos(int cidx, const std::string& op, amrex::Real check_val);

protected:
    /*void prepare_netcdf_file() override {}
    void process_output() override
    {
        // Test buffer populate for GPU runs
        std::vector<double> buf(num_gridpoints());
        //sampling_container().populate_buffer(buf);

    }*/
    const amrex::Real tol = 1e-8;
};

int FreeSurfaceImpl::check_output(const std::string& op, amrex::Real check_val)
{
    // Get number of points and output array
    auto npts_tot = num_gridpoints();
    auto out = heights();
    // Loop through grid points and check output
    int icheck = 0;
    for (int n = 0; n < npts_tot; ++n) {
        if (op == "=") {
            EXPECT_EQ(out[n], check_val);
            ++icheck;
        } else {
            if (op == "<") {
                EXPECT_LT(out[n], check_val);
                ++icheck;
            } else {
                if (op == "~") {
                    EXPECT_NEAR(out[n], check_val, tol);
                    ++icheck;
                }
            }
        }
    }
    return icheck;
}

int FreeSurfaceImpl::check_pos(
    const int cidx, const std::string& op, amrex::Real check_val)
{
    // Get number of points and position array
    auto npts_tot = num_gridpoints();
    auto locs = locations();
    // Loop through grid points and check output
    int icheck = 0;
    for (int n = 0; n < npts_tot; ++n) {
        if (op == "=") {
            EXPECT_EQ(locs[n][cidx], check_val);
            ++icheck;
        } else {
            if (op == "<") {
                EXPECT_LT(locs[n][cidx], check_val);
                ++icheck;
            } else {
                if (op == "~") {
                    EXPECT_NEAR(locs[n][cidx], check_val, tol);
                    ++icheck;
                }
            }
        }
    }
    return icheck;
}

} // namespace

class FreeSurfaceTest : public MeshTest
{
protected:
    void populate_parameters() override
    {
        MeshTest::populate_parameters();

        {
            amrex::ParmParse pp("amr");
            amrex::Vector<int> ncell{{32, 32, 64}};
            pp.add("max_level", 0);
            pp.addarr("n_cell", ncell);
        }
        {
            amrex::ParmParse pp("geometry");

            pp.addarr("prob_lo", problo);
            pp.addarr("prob_hi", probhi);
            pp.addarr("is_periodic", amrex::Vector<int>{{1, 1, 0}});
        }
    }
    void setup_grid0D()
    {
        amrex::ParmParse pp("freesurface");
        pp.add("output_frequency", 1);
        pp.addarr("num_points", amrex::Vector<int>{1, 1});
        pp.addarr("start", pt_coord);
        pp.addarr("end", pt_coord);
    }
    void setup_grid2D()
    {
        amrex::ParmParse pp("freesurface");
        pp.add("output_frequency", 1);
        pp.addarr("num_points", amrex::Vector<int>{npts, npts});
        pp.addarr("start", pl_start);
        pp.addarr("end", pl_end);
    }
    // Parameters to reuse
    const amrex::Real water_level0 = 64.0, water_level1 = 31.5;
    const amrex::Vector<amrex::Real> problo{{0.0, 0.0, -4.0}};
    const amrex::Vector<amrex::Real> probhi{{128.0, 128.0, 124.0}};
    const amrex::Vector<amrex::Real> pt_coord{{63.0, 65.0, 0.0}};
    const amrex::Vector<amrex::Real> pl_start{{0.0, 0.0, 0.0}};
    const amrex::Vector<amrex::Real> pl_end{{128.0, 128.0, 0.0}};
    const int npts = 3;
};

TEST_F(FreeSurfaceTest, point)
{
    initialize_mesh();
    auto& repo = sim().repo();
    auto& vof = repo.declare_field("vof", 1, 2);
    setup_grid0D();

    amrex::Real liwl = init_vof(vof, water_level0);
    auto& m_sim = sim();
    FreeSurfaceImpl tool(m_sim, "freesurface");
    tool.initialize();
    tool.post_advance_work();

    // Check number of points
    auto ngp = tool.num_gridpoints();
    EXPECT_EQ(ngp, 1);
    // Check location after being read
    int npos = tool.check_pos(0, "=", pt_coord[0]);
    ASSERT_EQ(npos, 1);
    npos = tool.check_pos(1, "=", pt_coord[1]);
    ASSERT_EQ(npos, 1);
    // Check output value
    int nout = tool.check_output("~", liwl);
    ASSERT_EQ(nout, 1);
}

TEST_F(FreeSurfaceTest, plane)
{
    initialize_mesh();
    auto& repo = sim().repo();
    auto& vof = repo.declare_field("vof", 1, 2);
    setup_grid2D();

    amrex::Real liwl = init_vof(vof, water_level1);
    auto& m_sim = sim();
    FreeSurfaceImpl tool(m_sim, "freesurface");
    tool.initialize();
    tool.post_advance_work();

    // Check number of points
    auto ngp = tool.num_gridpoints();
    EXPECT_EQ(ngp, npts * npts);
    // Check output value
    int nout = tool.check_output("~", liwl);
    ASSERT_EQ(nout, npts * npts);
}

} // namespace amr_wind_tests
