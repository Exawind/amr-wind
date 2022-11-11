#include "aw_test_utils/MeshTest.H"
#include "amr-wind/utilities/sampling/FreeSurface.H"
#include "amr-wind/utilities/tagging/FieldRefinement.H"

namespace amr_wind_tests {

namespace {

void init_vof(amr_wind::Field& vof_fld, amrex::Real water_level)
{
    const auto& mesh = vof_fld.repo().mesh();
    const int nlevels = vof_fld.repo().num_active_levels();

    // Since VOF is cell centered
    amrex::Real offset = 0.5;
    // VOF has only one component
    int d = 0;

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
    }
}

void init_vof_multival(
    amr_wind::Field& vof_fld, amrex::Real wl2, amrex::Real wl1, amrex::Real wl0)
{
    const auto& mesh = vof_fld.repo().mesh();
    const int nlevels = vof_fld.repo().num_active_levels();

    // This function initializes a vof field divided by a liquid-gas interface
    // at wl0, above which is another layer of liquid, introducing interfaces
    // at wl1 and wl2.

    // Since VOF is cell centered
    amrex::Real offset = 0.5;
    // VOF has only one component
    int d = 0;

    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& dx = mesh.Geom(lev).CellSizeArray();
        const auto& problo = mesh.Geom(lev).ProbLoArray();

        for (amrex::MFIter mfi(vof_fld(lev)); mfi.isValid(); ++mfi) {
            auto bx = mfi.growntilebox();
            const auto& farr = vof_fld(lev).array(mfi);

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                const amrex::Real z = problo[2] + (k + offset) * dx[2];
                amrex::Real local_vof;
                // Above wl1
                if (z - offset * dx[2] > wl1) {
                    local_vof = std::min(
                        1.0,
                        std::max(0.0, (wl2 - (z - offset * dx[2])) / dx[2]));
                } else {
                    // Above wl0
                    if (z - offset * dx[2] > wl0) {
                        local_vof = std::min(
                            1.0,
                            std::max(
                                0.0, ((z + offset * dx[2]) - wl1) / dx[2]));
                    } else {
                        // Bottom portion
                        local_vof = std::min(
                            1.0,
                            std::max(
                                0.0, (wl0 - (z - offset * dx[2])) / dx[2]));
                    }
                }
                farr(i, j, k, d) = local_vof;
            });
        }
    }
}

void init_vof_slope(
    amr_wind::Field& vof_fld,
    amrex::Real water_level,
    amrex::Real slope,
    amrex::Real domain_length)
{
    const auto& mesh = vof_fld.repo().mesh();
    const int nlevels = vof_fld.repo().num_active_levels();

    // The interface has a slope in x and y and intersects the center of the
    // domain at the designated water level Since VOF is cell centered
    amrex::Real offset = 0.5;
    // VOF has only one component
    int d = 0;

    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& dx = mesh.Geom(lev).CellSizeArray();
        const auto& problo = mesh.Geom(lev).ProbLoArray();

        for (amrex::MFIter mfi(vof_fld(lev)); mfi.isValid(); ++mfi) {
            auto bx = mfi.growntilebox();
            const auto& farr = vof_fld(lev).array(mfi);

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                const amrex::Real x = problo[0] + (i + offset) * dx[0];
                const amrex::Real y = problo[1] + (j + offset) * dx[1];
                const amrex::Real z = problo[2] + (k + offset) * dx[2];

                // Find height of interface at current x, y
                amrex::Real local_ht = water_level +
                                       slope * (x - 0.5 * domain_length) +
                                       slope * (y - 0.5 * domain_length);

                amrex::Real local_vof = std::min(
                    1.0,
                    std::max(0.0, (local_ht - (z - offset * dx[2])) / dx[2]));
                farr(i, j, k, d) = local_vof;
            });
        }
    }
}

//! Custom mesh class to be able to refine like a simulation would
//  - combination of AmrTestMesh and incflo classes
//  - with ability to initialize the refiner and regrid
class FSRefineMesh : public AmrTestMesh
{
public:
    FSRefineMesh()
        : m_sim(*this)
        , m_repo(m_sim.repo())
        , m_mesh_refiner(new amr_wind::RefineCriteriaManager(m_sim))
    {}
    amr_wind::FieldRepo& field_repo() { return m_repo; }
    amr_wind::CFDSim& sim() { return m_sim; }
    void init_refiner() { m_mesh_refiner->initialize(); }
    void remesh() { regrid(0, 0.0); }

protected:
    void MakeNewLevelFromScratch(
        int lev,
        amrex::Real time,
        const amrex::BoxArray& ba,
        const amrex::DistributionMapping& dm) override
    {
        SetBoxArray(lev, ba);
        SetDistributionMap(lev, dm);

        m_repo.make_new_level_from_scratch(lev, time, ba, dm);
    }

    void MakeNewLevelFromCoarse(
        int lev,
        amrex::Real time,
        const amrex::BoxArray& ba,
        const amrex::DistributionMapping& dm) override
    {
        SetBoxArray(lev, ba);
        SetDistributionMap(lev, dm);

        m_repo.make_new_level_from_coarse(lev, time, ba, dm);
    }

    void RemakeLevel(
        int lev,
        amrex::Real time,
        const amrex::BoxArray& ba,
        const amrex::DistributionMapping& dm) override
    {
        SetBoxArray(lev, ba);
        SetDistributionMap(lev, dm);

        m_repo.remake_level(lev, time, ba, dm);
    }

    void ClearLevel(int lev) override { m_repo.clear_level(lev); }

    void ErrorEst(
        int lev, amrex::TagBoxArray& tags, amrex::Real time, int ngrow) override
    {
        m_mesh_refiner->tag_cells(lev, tags, time, ngrow);
    }

private:
    amr_wind::CFDSim m_sim;
    amr_wind::FieldRepo& m_repo;
    std::unique_ptr<amr_wind::RefineCriteriaManager> m_mesh_refiner;
};

class FreeSurfaceImpl : public amr_wind::free_surface::FreeSurface
{
public:
    FreeSurfaceImpl(amr_wind::CFDSim& sim, const std::string& label)
        : amr_wind::free_surface::FreeSurface(sim, label)
    {}
    int check_output(const std::string& op, amrex::Real check_val);
    int check_output(int cidx, const std::string& op, amrex::Real check_val);
    int check_output_vec(
        const std::string& op, amrex::Vector<amrex::Real> check_val);
    int check_pos(int cidx, const std::string& op, amrex::Real check_val);

protected:
    // No file output during test
    void prepare_netcdf_file() override {}
    void process_output() override {}
    const amrex::Real tol = 1e-8;
};

int FreeSurfaceImpl::check_output(
    int cidx, const std::string& op, amrex::Real check_val)
{
    // Get number of points and output array
    auto npts_tot = num_gridpoints();
    auto out = heights();
    // Loop through grid points and check output
    int icheck = 0;
    for (int n = 0; n < npts_tot; ++n) {
        if (op == "=") {
            EXPECT_EQ(out[cidx * npts_tot + n], check_val);
            ++icheck;
        } else {
            if (op == "<") {
                EXPECT_LT(out[cidx * npts_tot + n], check_val);
                ++icheck;
            } else {
                if (op == "~") {
                    EXPECT_NEAR(out[cidx * npts_tot + n], check_val, tol);
                    ++icheck;
                }
            }
        }
    }
    return icheck;
}

int FreeSurfaceImpl::check_output(const std::string& op, amrex::Real check_val)
{
    return check_output(0, op, check_val);
}

int FreeSurfaceImpl::check_output_vec(
    const std::string& op, amrex::Vector<amrex::Real> check_val)
{
    // Get number of points and output array
    auto npts_tot = num_gridpoints();
    auto out = heights();
    // Loop through grid points and check output
    int icheck = 0;
    for (int n = 0; n < npts_tot; ++n) {
        if (op == "=") {
            EXPECT_EQ(out[n], check_val[n]);
            ++icheck;
        } else {
            if (op == "<") {
                EXPECT_LT(out[n], check_val[n]);
                ++icheck;
            } else {
                if (op == "~") {
                    EXPECT_NEAR(out[n], check_val[n], tol);
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
            pp.add("max_level", m_nlev);
            pp.addarr("n_cell", ncell);
        }
        {
            amrex::ParmParse pp("geometry");

            pp.addarr("prob_lo", problo);
            pp.addarr("prob_hi", probhi);
            pp.addarr("is_periodic", amrex::Vector<int>{{1, 1, 0}});
        }
    }
    void setup_grid0D(int ninst, std::string fsname)
    {
        amrex::ParmParse pp(fsname);
        pp.add("output_frequency", 1);
        pp.add("num_instances", ninst);
        pp.addarr("num_points", amrex::Vector<int>{1, 1});
        pp.addarr("start", pt_coord);
        pp.addarr("end", pt_coord);
    }
    void setup_grid0D(int ninst) { setup_grid0D(ninst, "freesurface"); }
    void setup_grid2D(int ninst)
    {
        amrex::ParmParse pp("freesurface");
        pp.add("output_frequency", 1);
        pp.add("num_instances", ninst);
        pp.addarr("num_points", amrex::Vector<int>{npts, npts});
        pp.addarr("start", pl_start);
        pp.addarr("end", pl_end);
    }
    void setup_grid2D_narrow()
    {
        amrex::ParmParse pp("freesurface");
        pp.add("output_frequency", 1);
        pp.add("num_instances", 1);
        pp.addarr("num_points", amrex::Vector<int>{npts, npts});
        pp.addarr("start", plnarrow_s);
        pp.addarr("end", plnarrow_e);
    }
    void setup_fieldrefinement()
    {
        amrex::ParmParse pp("tagging");
        pp.add("labels", (std::string) "t1");
        amrex::ParmParse ppt1("tagging.t1");
        ppt1.add("type", (std::string) "FieldRefinement");
        ppt1.add("field_name", fname);
        ppt1.addarr("field_error", amrex::Vector<amrex::Real>{fref_val});
    }
    // Parameters to reuse
    const amrex::Real water_level0 = 64.0;
    const amrex::Real water_level1 = 31.5;
    const amrex::Real water_level2 = 65.0;
    const amrex::Vector<amrex::Real> problo{{0.0, 0.0, -4.0}};
    const amrex::Vector<amrex::Real> probhi{{128.0, 128.0, 124.0}};
    const amrex::Vector<amrex::Real> pt_coord{{63.0, 65.0, 0.0}};
    const amrex::Vector<amrex::Real> pl_start{{0.0, 0.0, 0.0}};
    const amrex::Vector<amrex::Real> pl_end{{128.0, 128.0, 0.0}};
    const amrex::Vector<amrex::Real> plnarrow_s{{63.0, 63.0, 0.0}};
    const amrex::Vector<amrex::Real> plnarrow_e{{65.0, 65.0, 0.0}};
    const int npts = 3;
    const amrex::Real fref_val = 0.5;
    const std::string fname = "flag";
    int m_nlev = 0;
};

TEST_F(FreeSurfaceTest, point)
{
    initialize_mesh();
    auto& repo = sim().repo();
    auto& vof = repo.declare_field("vof", 1, 2);
    setup_grid0D(1);

    init_vof(vof, water_level0);
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
    int nout = tool.check_output("~", water_level0);
    ASSERT_EQ(nout, 1);
}

TEST_F(FreeSurfaceTest, plane)
{
    initialize_mesh();
    auto& repo = sim().repo();
    auto& vof = repo.declare_field("vof", 1, 2);
    setup_grid2D(1);

    init_vof(vof, water_level1);
    auto& m_sim = sim();
    FreeSurfaceImpl tool(m_sim, "freesurface");
    tool.initialize();
    tool.post_advance_work();

    // Check number of points
    auto ngp = tool.num_gridpoints();
    EXPECT_EQ(ngp, npts * npts);
    // Check output value
    int nout = tool.check_output("~", water_level1);
    ASSERT_EQ(nout, npts * npts);
}

TEST_F(FreeSurfaceTest, multivalued)
{
    initialize_mesh();
    auto& repo = sim().repo();
    auto& vof = repo.declare_field("vof", 1, 2);
    setup_grid2D(3);

    // Parameters of water level
    amrex::Real wl0 = probhi[2] * 2 / 3;
    amrex::Real wl1 = probhi[2] / 2;
    amrex::Real wl2 = probhi[2] / 5;

    init_vof_multival(vof, wl0, wl1, wl2);
    auto& m_sim = sim();
    FreeSurfaceImpl tool(m_sim, "freesurface");
    tool.initialize();
    tool.post_advance_work();

    // Check number of outputs
    auto heights = tool.heights();
    EXPECT_EQ(heights.size(), 3 * npts * npts);

    // Check output values
    int nout = tool.check_output(0, "~", wl0);
    ASSERT_EQ(nout, npts * npts);
    nout = tool.check_output(1, "~", wl1);
    ASSERT_EQ(nout, npts * npts);
    nout = tool.check_output(2, "~", wl2);
    ASSERT_EQ(nout, npts * npts);
}

TEST_F(FreeSurfaceTest, sloped)
{
    initialize_mesh();
    auto& repo = sim().repo();
    auto& vof = repo.declare_field("vof", 1, 2);
    setup_grid2D_narrow();

    amrex::Real slope = 0.125;
    amrex::Real domain_l = probhi[0];
    init_vof_slope(vof, water_level2, slope, domain_l);
    auto& m_sim = sim();
    FreeSurfaceImpl tool(m_sim, "freesurface");
    tool.initialize();
    tool.post_advance_work();

    // Calculate expected output values
    amrex::Vector<amrex::Real> out_vec;
    out_vec.resize(npts * npts, 0.0);
    // Step in x, then y
    out_vec[0] = water_level2 + slope * (-1.0 - 1.0);
    out_vec[1] = water_level2 + slope * (+0.0 - 1.0);
    out_vec[2] = water_level2 + slope * (+1.0 - 1.0);
    out_vec[3] = water_level2 + slope * (-1.0 + 0.0);
    out_vec[4] = water_level2 + slope * (+0.0 + 0.0);
    out_vec[5] = water_level2 + slope * (+1.0 + 0.0);
    out_vec[6] = water_level2 + slope * (-1.0 + 1.0);
    out_vec[7] = water_level2 + slope * (+0.0 + 1.0);
    out_vec[8] = water_level2 + slope * (+1.0 + 1.0);
    // Check output value
    int nout = tool.check_output_vec("~", out_vec);
    ASSERT_EQ(nout, npts * npts);
}

TEST_F(FreeSurfaceTest, multisampler)
{
    initialize_mesh();
    auto& repo = sim().repo();
    auto& vof = repo.declare_field("vof", 1, 2);

    // Set up parameters for one sampler
    setup_grid0D(1, "freesurface0");

    // Set up parameters for another sampler
    setup_grid2D_narrow();

    // Initialize VOF distribution and access sim
    init_vof(vof, water_level1);
    auto& m_sim = sim();

    // Initialize first sampler
    FreeSurfaceImpl tool1(m_sim, "freesurface0");
    tool1.initialize();

    // Initialize second sampler
    FreeSurfaceImpl tool2(m_sim, "freesurface");
    tool2.initialize();
}

TEST_F(FreeSurfaceTest, regrid)
{
    // Allow 2 levels
    m_nlev = 1;
    // Set up parameters for domain
    populate_parameters();
    // Set up parameters for refinement
    setup_fieldrefinement();
    // Set up parameters for sampler
    setup_grid2D(1);
    // Create mesh and initialize
    reset_prob_domain();
    auto rmesh = FSRefineMesh();
    rmesh.initialize_mesh(0.0);

    // Repo and fields
    auto& repo = rmesh.field_repo();
    auto& vof = repo.declare_field("vof", 1, 2);
    auto& flag = repo.declare_field(fname, 1, 2);

    // Set up scalar for determining refinement - all fine level
    flag.setVal(2.0 * fref_val);

    // Initialize mesh refiner and remesh
    rmesh.init_refiner();
    rmesh.remesh();

    // Initialize VOF distribution and access sim
    init_vof(vof, water_level1);
    auto& rsim = rmesh.sim();

    // Initialize sampler and check result on initial mesh
    FreeSurfaceImpl tool(rsim, "freesurface");
    tool.initialize();
    tool.post_advance_work();
    tool.check_output("~", water_level1);

    // Change scalar for determining refinement - no fine level
    flag.setVal(0.0);

    // Regrid, update fields, and do post_regrid_actions for tool
    rmesh.remesh();
    tool.post_regrid_actions();

    // Check that result is unchanged on new mesh
    tool.post_advance_work();
    tool.check_output("~", water_level1);
}

} // namespace amr_wind_tests
