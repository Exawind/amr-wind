#include "aw_test_utils/MeshTest.H"
#include "amr-wind/utilities/IOManager.H"
#include "amr-wind/utilities/sampling/FieldNorms.H"
#include "amr-wind/utilities/tagging/FieldRefinement.H"

namespace amr_wind_tests {

namespace {

void init_velocity(
    amr_wind::Field& vel_fld,
    amrex::Real u,
    amrex::Real v,
    amrex::Real w,
    amrex::Real lev0_fac)
{
    const int nlevels = vel_fld.repo().num_active_levels();
    const amrex::GpuArray<amrex::Real, 3> vels = {u, v, w};
    for (int lev = 0; lev < nlevels; ++lev) {

        // Apply a factor to the base level values
        const amrex::Real fac = (lev == 0) ? lev0_fac : 1.0;
        const auto& farrs = vel_fld(lev).arrays();

        amrex::ParallelFor(
            vel_fld(lev), vel_fld.num_grow(), vel_fld.num_comp(),
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k, int n) noexcept {
                // Mix positive and negative as check on L2 norm
                farrs[nbx](i, j, k, n) =
                    fac * (i % 2 == 0 ? vels[n] : -vels[n]);
            });
    }
    amrex::Gpu::synchronize();
}

//! Custom mesh class to be able to refine like a simulation would
//  - combination of AmrTestMesh and incflo classes
//  - with ability to initialize the refiner and regrid
class FNRefinemesh : public AmrTestMesh
{
public:
    FNRefinemesh() : m_mesh_refiner(new amr_wind::RefineCriteriaManager(m_sim))
    {}
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
    std::unique_ptr<amr_wind::RefineCriteriaManager> m_mesh_refiner;
};

class FieldNormsImpl : public amr_wind::field_norms::FieldNorms
{
public:
    FieldNormsImpl(amr_wind::CFDSim& sim, const std::string& label)
        : amr_wind::field_norms::FieldNorms(sim, label)
    {}
    void check_output(
        amrex::Real check_val0, amrex::Real check_val1, amrex::Real check_val2);

protected:
    // No file output during test
    void prepare_ascii_file() override {}
    void write_ascii() override {}
    const amrex::Real m_tol = 1e-8;
};

void FieldNormsImpl::check_output(
    amrex::Real check_val0, amrex::Real check_val1, amrex::Real check_val2)
{
    // Get variable names and check
    EXPECT_EQ(var_names()[0], (std::string) "velocityx");
    EXPECT_EQ(var_names()[1], (std::string) "velocityy");
    EXPECT_EQ(var_names()[2], (std::string) "velocityz");
    // Loop through norm values and check them
    const amrex::Real tiny = std::numeric_limits<amrex::Real>::epsilon();
    EXPECT_NEAR(field_norms()[0], check_val0, 500 * tiny);
    EXPECT_NEAR(field_norms()[1], check_val1, 500 * tiny);
    EXPECT_NEAR(field_norms()[2], check_val2, 500 * tiny);
}

} // namespace

class FieldNormsTest : public MeshTest
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

            pp.addarr("prob_lo", m_problo);
            pp.addarr("prob_hi", m_probhi);
            pp.addarr("is_periodic", amrex::Vector<int>{{1, 1, 0}});
        }
    }
    static void setup_fnorm(bool levelmask_flag)
    {
        amrex::ParmParse pp("fieldnorm");
        pp.add("output_frequency", 1);
        if (!levelmask_flag) {
            pp.add("mask_redundant_grids", false);
        }
    }
    void setup_fieldrefinement()
    {
        amrex::ParmParse pp("tagging");
        pp.add("labels", (std::string) "t1");
        amrex::ParmParse ppt1("tagging.t1");
        ppt1.add("type", (std::string) "FieldRefinement");
        ppt1.add("field_name", m_fname);
        ppt1.addarr("field_error", amrex::Vector<amrex::Real>{m_fref_val});
    }
    void setup_intfieldrefinement()
    {
        amrex::ParmParse pp("tagging");
        pp.add("labels", (std::string) "t1");
        amrex::ParmParse ppt1("tagging.t1");
        ppt1.add("type", (std::string) "FieldRefinement");
        ppt1.add("field_name", m_ifname);
        ppt1.addarr("field_error", amrex::Vector<int>{m_ifref_val});
    }
    static void setup_tagging_box()
    {
        amrex::ParmParse ppt1("tagging.t1");
        amrex::Vector<amrex::Real> blo = {3.0, 3.0, 3.0};
        amrex::Vector<amrex::Real> bhi = {30.0, 30.0, 30.0};
        ppt1.addarr("box_lo", blo);
        ppt1.addarr("box_hi", bhi);
    }
    // Parameters to reuse
    const amrex::Vector<amrex::Real> m_problo{{0.0, 0.0, -4.0}};
    const amrex::Vector<amrex::Real> m_probhi{{128.0, 128.0, 124.0}};
    const amrex::Real m_fref_val = 0.5;
    const std::string m_fname = "flag";
    const int m_ifref_val = 1;
    const std::string m_ifname = "iflag";
    const int m_nlev = 1;
    const int m_nx = 32;
    const int m_ny = 32;
    const int m_nz = 64;
    // Number of cells on base level
    const int m_ncell0 = m_nx * m_ny * m_nz;
    // Cell volume on base level
    const amrex::Real m_cv0 = (128.0 / m_nx) * (128.0 / m_ny) * (128.0 / m_nz);
    // Total domain volume
    const amrex::Real m_dv = 128. * 128. * 128.;
    // Velocity component values
    const amrex::Real m_u = 0.9;
    const amrex::Real m_v = 1.3;
    const amrex::Real m_w = 2.7;
};

TEST_F(FieldNormsTest, levelmask_on)
{
    bool levelmask = true;
    // Set up parameters for domain
    populate_parameters();
    // Set up parameters for refinement
    setup_fieldrefinement();
    // Set up parameters for sampler
    setup_fnorm(levelmask);
    // Create mesh and initialize
    reset_prob_domain();
    auto rmesh = FNRefinemesh();
    rmesh.initialize_mesh(0.0);

    // Repo and fields
    auto& repo = rmesh.field_repo();
    auto& velocity = repo.declare_field("velocity", 3, 2);
    auto& flag = repo.declare_field(m_fname, 1, 2);

    // Set up scalar for determining refinement - all fine level
    flag.setVal(2.0 * m_fref_val);

    // Initialize mesh refiner and remesh
    rmesh.init_refiner();
    rmesh.remesh();

    // Initialize velocity distribution and access sim
    const amrex::Real lev0_fac = 1.5;
    init_velocity(velocity, m_u, m_v, m_w, lev0_fac);
    auto& rsim = rmesh.sim();

    // Initialize IOManager because FieldNorms relies on it
    auto& io_mgr = rsim.io_manager();
    // Set up velocity as an output (plot) variable
    io_mgr.register_output_var("velocity");
    io_mgr.initialize_io();

    // Initialize sampler and check result on initial mesh
    FieldNormsImpl tool(rsim, "fieldnorm");
    tool.initialize();
    tool.post_advance_work();
    // Only highest level will be counted
    // sqrt(number of cells * cell volume * cell value * cell value)
    amrex::Real unorm =
        std::sqrt((m_ncell0 * 8.) * (m_cv0 / 8.) * m_u * m_u / m_dv);
    amrex::Real vnorm =
        std::sqrt((m_ncell0 * 8.) * (m_cv0 / 8.) * m_v * m_v / m_dv);
    amrex::Real wnorm =
        std::sqrt((m_ncell0 * 8.) * (m_cv0 / 8.) * m_w * m_w / m_dv);
    tool.check_output(unorm, vnorm, wnorm);

    // Change scalar for determining refinement - no fine level
    flag.setVal(0.0);

    // Regrid
    rmesh.remesh();

    // Check result on new mesh
    tool.post_advance_work();
    // Only base level exists, which includes factor
    unorm =
        std::sqrt(m_ncell0 * m_cv0 * lev0_fac * lev0_fac * m_u * m_u / m_dv);
    vnorm =
        std::sqrt(m_ncell0 * m_cv0 * lev0_fac * lev0_fac * m_v * m_v / m_dv);
    wnorm =
        std::sqrt(m_ncell0 * m_cv0 * lev0_fac * lev0_fac * m_w * m_w / m_dv);
    tool.check_output(unorm, vnorm, wnorm);
}

TEST_F(FieldNormsTest, levelmask_on_with_box)
{
    bool levelmask = true;
    // Set up parameters for domain
    populate_parameters();
    // Set up parameters for refinement
    setup_fieldrefinement();
    setup_tagging_box();
    // Set up parameters for sampler
    setup_fnorm(levelmask);
    // Create mesh and initialize
    reset_prob_domain();
    auto rmesh = FNRefinemesh();
    rmesh.initialize_mesh(0.0);

    // Repo and fields
    auto& repo = rmesh.field_repo();
    auto& velocity = repo.declare_field("velocity", 3, 2);
    auto& flag = repo.declare_field(m_fname, 1, 2);

    // Set up scalar for determining refinement - all fine level
    flag.setVal(2.0 * m_fref_val);

    // Initialize mesh refiner and remesh
    rmesh.init_refiner();
    rmesh.remesh();

    // Initialize velocity distribution and access sim
    const amrex::Real lev0_fac = 1.5;
    init_velocity(velocity, m_u, m_v, m_w, lev0_fac);
    auto& rsim = rmesh.sim();

    // Initialize IOManager because FieldNorms relies on it
    auto& io_mgr = rsim.io_manager();
    // Set up velocity as an output (plot) variable
    io_mgr.register_output_var("velocity");
    io_mgr.initialize_io();

    // Initialize sampler and check result on initial mesh
    FieldNormsImpl tool(rsim, "fieldnorm");
    tool.initialize();
    tool.post_advance_work();
    tool.check_output(
        1.342655804506534, 1.9393917176204616, 4.0279674135195087);

    // Change scalar for determining refinement - no fine level
    flag.setVal(0.0);

    // Regrid
    rmesh.remesh();

    // Check result on new mesh
    tool.post_advance_work();
    tool.check_output(
        1.3500000000000343, 1.9499999999999589, 4.0500000000000043);
}

// Not sure why this option would be used, but it is available
TEST_F(FieldNormsTest, levelmask_off)
{
    bool levelmask = false;
    // Set up parameters for domain
    populate_parameters();
    // Set up parameters for refinement
    setup_fieldrefinement();
    // Set up parameters for sampler
    setup_fnorm(levelmask);
    // Create mesh and initialize
    reset_prob_domain();
    auto rmesh = FNRefinemesh();
    rmesh.initialize_mesh(0.0);

    // Repo and fields
    auto& repo = rmesh.field_repo();
    auto& velocity = repo.declare_field("velocity", 3, 2);
    auto& flag = repo.declare_field(m_fname, 1, 2);

    // Set up scalar for determining refinement - all fine level
    flag.setVal(2.0 * m_fref_val);

    // Initialize mesh refiner and remesh
    rmesh.init_refiner();
    rmesh.remesh();

    // Initialize velocity distribution and access sim
    const amrex::Real lev0_fac = 1.5;
    init_velocity(velocity, m_u, m_v, m_w, lev0_fac);
    auto& rsim = rmesh.sim();

    // Initialize IOManager because FieldNorms relies on it
    auto& io_mgr = rsim.io_manager();
    // Set up velocity as an output (plot) variable
    io_mgr.register_output_var("velocity");
    io_mgr.initialize_io();

    // Initialize sampler and check result on initial mesh
    FieldNormsImpl tool(rsim, "fieldnorm");
    tool.initialize();
    tool.post_advance_work();
    // Both levels will be counted
    amrex::Real unorm = std::sqrt(
        ((m_ncell0 * 8.) * (m_cv0 / 8.) * m_u * m_u +
         m_ncell0 * m_cv0 * lev0_fac * lev0_fac * m_u * m_u) /
        m_dv);
    amrex::Real vnorm = std::sqrt(
        ((m_ncell0 * 8.) * (m_cv0 / 8.) * m_v * m_v +
         m_ncell0 * m_cv0 * lev0_fac * lev0_fac * m_v * m_v) /
        m_dv);
    amrex::Real wnorm = std::sqrt(
        ((m_ncell0 * 8.) * (m_cv0 / 8.) * m_w * m_w +
         m_ncell0 * m_cv0 * lev0_fac * lev0_fac * m_w * m_w) /
        m_dv);
    tool.check_output(unorm, vnorm, wnorm);

    // Change scalar for determining refinement - no fine level
    flag.setVal(0.0);

    // Regrid
    rmesh.remesh();

    // Check result on new mesh
    tool.post_advance_work();
    // Only base level exists, which includes factor
    unorm =
        std::sqrt(m_ncell0 * m_cv0 * lev0_fac * lev0_fac * m_u * m_u / m_dv);
    vnorm =
        std::sqrt(m_ncell0 * m_cv0 * lev0_fac * lev0_fac * m_v * m_v / m_dv);
    wnorm =
        std::sqrt(m_ncell0 * m_cv0 * lev0_fac * lev0_fac * m_w * m_w / m_dv);
    tool.check_output(unorm, vnorm, wnorm);
}

TEST_F(FieldNormsTest, levelmask_on_int_refinement)
{
    bool levelmask = true;
    // Set up parameters for domain
    populate_parameters();
    // Set up parameters for refinement
    setup_intfieldrefinement();
    // Set up parameters for sampler
    setup_fnorm(levelmask);
    // Create mesh and initialize
    reset_prob_domain();
    auto rmesh = FNRefinemesh();
    rmesh.initialize_mesh(0.0);

    // Repo and fields
    auto& repo = rmesh.field_repo();
    auto& velocity = repo.declare_field("velocity", 3, 2);
    auto& flag = repo.declare_int_field(m_ifname, 1, 2);

    // Set up scalar for determining refinement - all fine level
    flag.setVal(2 * m_ifref_val);

    // Initialize mesh refiner and remesh
    rmesh.init_refiner();
    rmesh.remesh();

    // Initialize velocity distribution and access sim
    const amrex::Real lev0_fac = 1.5;
    init_velocity(velocity, m_u, m_v, m_w, lev0_fac);
    auto& rsim = rmesh.sim();

    // Initialize IOManager because FieldNorms relies on it
    auto& io_mgr = rsim.io_manager();
    // Set up velocity as an output (plot) variable
    io_mgr.register_output_var("velocity");
    io_mgr.initialize_io();

    // Initialize sampler and check result on initial mesh
    FieldNormsImpl tool(rsim, "fieldnorm");
    tool.initialize();
    tool.post_advance_work();
    // Only highest level will be counted
    // sqrt(number of cells * cell volume * cell value * cell value)
    amrex::Real unorm =
        std::sqrt((m_ncell0 * 8.) * (m_cv0 / 8.) * m_u * m_u / m_dv);
    amrex::Real vnorm =
        std::sqrt((m_ncell0 * 8.) * (m_cv0 / 8.) * m_v * m_v / m_dv);
    amrex::Real wnorm =
        std::sqrt((m_ncell0 * 8.) * (m_cv0 / 8.) * m_w * m_w / m_dv);
    tool.check_output(unorm, vnorm, wnorm);

    // Change scalar for determining refinement - no fine level
    flag.setVal(0.0);

    // Regrid
    rmesh.remesh();

    // Check result on new mesh
    tool.post_advance_work();
    // Only base level exists, which includes factor
    unorm =
        std::sqrt(m_ncell0 * m_cv0 * lev0_fac * lev0_fac * m_u * m_u / m_dv);
    vnorm =
        std::sqrt(m_ncell0 * m_cv0 * lev0_fac * lev0_fac * m_v * m_v / m_dv);
    wnorm =
        std::sqrt(m_ncell0 * m_cv0 * lev0_fac * lev0_fac * m_w * m_w / m_dv);
    tool.check_output(unorm, vnorm, wnorm);
}

} // namespace amr_wind_tests
