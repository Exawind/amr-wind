
#include "aw_test_utils/MeshTest.H"

#include "amr-wind/utilities/sampling/Sampling.H"
#include "amr-wind/utilities/sampling/SamplingContainer.H"
#include "amr-wind/utilities/sampling/ProbeSampler.H"
#include "amr-wind/utilities/sampling/PlaneSampler.H"
#include "amr-wind/utilities/sampling/VolumeSampler.H"
#include "amr-wind/utilities/sampling/DTUSpinnerSampler.H"
#include "amr-wind/utilities/sampling/RadarSampler.H"
#include "amr-wind/utilities/sampling/SamplingUtils.H"
#include "AMReX_Vector.H"
#include "amr-wind/core/vs/vector_space.H"
#include "amr-wind/utilities/tensor_ops.H"

namespace amr_wind_tests {

namespace {

void init_field(amr_wind::Field& fld)
{
    const auto& mesh = fld.repo().mesh();
    const int nlevels = fld.repo().num_active_levels();
    const int ncomp = fld.num_comp();

    amrex::Real offset = 0.0;
    if (fld.field_location() == amr_wind::FieldLoc::CELL) {
        offset = 0.5;
    }

    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& dx = mesh.Geom(lev).CellSizeArray();
        const auto& problo = mesh.Geom(lev).ProbLoArray();
        const auto& farrs = fld(lev).arrays();

        amrex::ParallelFor(
            fld(lev), fld.num_grow(), ncomp,
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k, int n) noexcept {
                const amrex::Real x = problo[0] + (i + offset) * dx[0];
                const amrex::Real y = problo[1] + (j + offset) * dx[1];
                const amrex::Real z = problo[2] + (k + offset) * dx[2];
                farrs[nbx](i, j, k, n) = x + y + z;
            });
    }
    amrex::Gpu::synchronize();
}

void init_int_field(amr_wind::IntField& fld)
{
    const int nlevels = fld.repo().num_active_levels();
    const int ncomp = fld.num_comp();

    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& farrs = fld(lev).arrays();

        amrex::ParallelFor(
            fld(lev), fld.num_grow(), ncomp,
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k, int n) noexcept {
                farrs[nbx](i, j, k, n) = i + j + k;
            });
    }
    amrex::Gpu::synchronize();
}

void write_probe_sampler_file(const std::string& fname)
{
    std::ofstream os(fname);
    // Total number of points
    os << "3\n";
    // Coordinates
    os << "0.0\t0.0\t0.0\n";
    os << "60.0\t2.0\t3.0\n";
    os << "100.0\t8.0\t5.0\n";
}

class SamplingImpl : public amr_wind::sampling::Sampling
{
public:
    SamplingImpl(amr_wind::CFDSim& sim, const std::string& label)
        : amr_wind::sampling::Sampling(sim, label)
    {}

    bool write_flag{false};

protected:
    void prepare_netcdf_file() override {}
    void process_output() override
    {
        // Test buffer populate for GPU runs
        std::vector<double> buf(
            num_total_particles() * var_names().size(), 0.0);
        sampling_container().populate_buffer(buf);

        write_flag = true;
    }
};

} // namespace

class SamplingTest : public MeshTest
{
protected:
    void populate_parameters() override
    {
        MeshTest::populate_parameters();

        {
            amrex::ParmParse pp("amr");
            amrex::Vector<int> ncell{{32, 32, 64}};
            pp.add("max_level", 0);
            pp.add("max_grid_size", 16);
            pp.addarr("n_cell", ncell);
        }
        {
            amrex::ParmParse pp("geometry");
            amrex::Vector<amrex::Real> problo{{0.0, 0.0, 0.0}};
            amrex::Vector<amrex::Real> probhi{{128.0, 128.0, 128.0}};

            pp.addarr("prob_lo", problo);
            pp.addarr("prob_hi", probhi);
        }
    }
};

namespace {

void test_scontainer_impl(
    amr_wind::sampling::SamplingContainer::ParticleVector& pvec, const int npts)
{
    using IIx = amr_wind::sampling::IIx;
    using PType = amr_wind::sampling::SamplingContainer::ParticleType;
    // Create a line probe
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> begin{{66.0, 66.0, 1.0}};
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> end{{66.0, 66.0, 127.0}};

    auto* pstruct = pvec.data();
    const int id_start = static_cast<int>(PType::NextID());
    const int proc_id = amrex::ParallelDescriptor::MyProc();

    amrex::ParallelFor(npts, [=] AMREX_GPU_DEVICE(const int ip) noexcept {
        auto& pp = pstruct[ip];
        pp.id() = id_start + ip;
        pp.cpu() = proc_id;
        pp.pos(0) = begin[0];
        pp.pos(1) = begin[1];

        const int dl = static_cast<int>(
            (end[2] - begin[2]) / static_cast<amrex::Real>(npts - 1));
        const amrex::Real z = (ip + 0.5) * dl;
        pp.pos(2) = z;

        pp.idata(IIx::uid) = static_cast<int>(pp.id());
        pp.idata(IIx::sid) = 0;
        pp.idata(IIx::nid) = ip;
    });
    amrex::Gpu::streamSynchronize();
}

} // namespace

TEST_F(SamplingTest, scontainer)
{
    if (amrex::ParallelDescriptor::NProcs() > 1) {
        GTEST_SKIP();
    }

    using ParIter = amr_wind::sampling::SamplingContainer::ParIterType;
    const int ncomp = 5;
    initialize_mesh();

    amr_wind::sampling::SamplingContainer sc(mesh());
    // Set up runtime components and prime the particle container
    sc.setup_container(ncomp);

    const int npts = 64;
    const int lev = 0;
    const int gid = 0;
    const int tid = 0;
    auto& ptile = sc.GetParticles(lev)[std::make_pair(gid, tid)];
    ptile.resize(npts);
    test_scontainer_impl(ptile.GetArrayOfStructs()(), npts);
    sc.Redistribute();

#ifndef AMREX_USE_GPU
    // TODO: Check why this causes errors in non-managed CUDA run
    ASSERT_EQ(npts, sc.NumberOfParticlesAtLevel(lev));
#endif
    int counter = 0;
    int total_particles = 0;
    for (ParIter pti(sc, lev); pti.isValid(); ++pti) {
        ++counter;
        total_particles += pti.numParticles();
    }
    ASSERT_EQ(total_particles, npts);
    ASSERT_EQ(counter, 4);
}

TEST_F(SamplingTest, sampling)
{
    initialize_mesh();
    auto& repo = sim().repo();
    auto& vel = repo.declare_field("velocity", 3, 2);
    auto& pres = repo.declare_nd_field("pressure", 1, 2);
    auto& rho = repo.declare_field("density", 1, 2);
    auto& idxsum = repo.declare_int_field("idxsum", 1, 2);
    init_field(vel);
    init_field(pres);
    init_field(rho);
    init_int_field(idxsum);

    {
        amrex::ParmParse pp("sampling");
        pp.add("output_frequency", 1);
        pp.addarr("labels", amrex::Vector<std::string>{"line1"});
        pp.addarr(
            "fields",
            amrex::Vector<std::string>{"density", "pressure", "velocity"});
        pp.addarr(
            "derived_fields", amrex::Vector<std::string>{"mag_vorticity"});
        pp.addarr("int_fields", amrex::Vector<std::string>{"idxsum"});
    }
    {
        amrex::ParmParse pp("sampling.line1");
        pp.add("type", std::string("LineSampler"));
        pp.add("num_points", 16);
        pp.addarr("start", amrex::Vector<amrex::Real>{66.0, 66.0, 1.0});
        pp.addarr("end", amrex::Vector<amrex::Real>{66.0, 66.0, 127.0});
    }

    SamplingImpl probes(sim(), "sampling");
    probes.initialize();
    probes.post_advance_work();

    EXPECT_TRUE(probes.write_flag);
}

TEST_F(SamplingTest, sampling_timing)
{
    initialize_mesh();
    auto& repo = sim().repo();
    auto& vel = repo.declare_field("velocity", 3, 2);
    auto& pres = repo.declare_nd_field("pressure", 1, 2);
    auto& rho = repo.declare_field("density", 1, 2);
    auto& idxsum = repo.declare_int_field("idxsum", 1, 2);
    init_field(vel);
    init_field(pres);
    init_field(rho);
    init_int_field(idxsum);

    {
        amrex::ParmParse pp("sampling");
        pp.add("output_frequency", 1);
        pp.add("output_delay", 1);
        pp.addarr("labels", amrex::Vector<std::string>{"line1"});
        pp.addarr(
            "fields",
            amrex::Vector<std::string>{"density", "pressure", "velocity"});
        pp.addarr(
            "derived_fields", amrex::Vector<std::string>{"mag_vorticity"});
        pp.addarr("int_fields", amrex::Vector<std::string>{"idxsum"});
    }
    {
        amrex::ParmParse pp("sampling.line1");
        pp.add("type", std::string("LineSampler"));
        pp.add("num_points", 16);
        pp.addarr("start", amrex::Vector<amrex::Real>{66.0, 66.0, 1.0});
        pp.addarr("end", amrex::Vector<amrex::Real>{66.0, 66.0, 127.0});
    }

    SamplingImpl probes(sim(), "sampling");
    probes.initialize();
    probes.post_advance_work();
    EXPECT_FALSE(probes.write_flag);
    sim().time().new_timestep();
    probes.post_advance_work();
    EXPECT_TRUE(probes.write_flag);
}

TEST_F(SamplingTest, probe_sampler)
{
    initialize_mesh();

    constexpr amrex::Real tol = 1.0e-12;
    std::string fname = "probes.txt";
    // Write file
    write_probe_sampler_file(fname);

    {
        amrex::ParmParse pp("cloud");
        pp.add("probe_location_file", fname);
        pp.addarr("offsets", amrex::Vector<double>{1.0, 2.5});
        pp.addarr("offset_vector", amrex::Vector<double>{0.2, 0.5, 1.0});
    }

    amr_wind::sampling::ProbeSampler cloud(sim());
    cloud.initialize("cloud");
    amr_wind::sampling::ProbeSampler::SampleLocType locs;
    cloud.sampling_locations(locs);

    ASSERT_EQ(locs.size(), 3 * 2);
    const amrex::Vector<amrex::Real> xprobe_golds{0.2, 60.2, 100.2,
                                                  0.5, 60.5, 100.5};
    const amrex::Vector<amrex::Real> yprobe_golds{0.5,  2.5,  8.5,
                                                  1.25, 3.25, 9.25};
    const amrex::Vector<amrex::Real> zprobe_golds{1.0, 4.0, 6.0, 2.5, 5.5, 7.5};
    for (int n = 0; n < locs.size(); ++n) {
        EXPECT_NEAR(locs[n][0], xprobe_golds[n], tol);
        EXPECT_NEAR(locs[n][1], yprobe_golds[n], tol);
        EXPECT_NEAR(locs[n][2], zprobe_golds[n], tol);
    }

    // Remove file
    const char* fname_char = fname.c_str();
    {
        std::ifstream f(fname_char);
        if (f.good()) {
            remove(fname_char);
        }
        // Check that file is removed
        std::ifstream ff(fname_char);
        EXPECT_FALSE(ff.good());
    }
}

TEST_F(SamplingTest, plane_sampler)
{
    initialize_mesh();

    {
        amrex::ParmParse pp("plane");
        pp.addarr("axis1", amrex::Vector<double>{0.0, 1.0, 0.0});
        pp.addarr("axis2", amrex::Vector<double>{0.0, 0.0, 1.0});
        pp.addarr("origin", amrex::Vector<double>{1.0, 1.0, 1.0});
        pp.addarr("num_points", amrex::Vector<int>{3, 3});
        pp.addarr("offsets", amrex::Vector<double>{2.0, 10.0});
        pp.addarr("offset_vector", amrex::Vector<double>{1.0, 0.0, 0.0});
    }

    amr_wind::sampling::PlaneSampler plane(sim());
    plane.initialize("plane");
    amr_wind::sampling::PlaneSampler::SampleLocType locs;
    plane.sampling_locations(locs);

    ASSERT_EQ(locs.size(), 3 * 3 * 2);
}

TEST_F(SamplingTest, volume_sampler)
{
    initialize_mesh();
    amrex::ParmParse pp("volume");
    pp.addarr("hi", amrex::Vector<double>{1.0, 1.0, 1.0});
    pp.addarr("lo", amrex::Vector<double>{0.0, 0.0, 0.0});
    pp.addarr("num_points", amrex::Vector<int>{3, 5, 5});

    amr_wind::sampling::VolumeSampler volume(sim());
    volume.initialize("volume");
    amr_wind::sampling::VolumeSampler::SampleLocType locs;
    volume.sampling_locations(locs);

    ASSERT_EQ(locs.size(), 3 * 5 * 5);
}

TEST_F(SamplingTest, spinner_sampler)
{
    initialize_mesh();
    amrex::ParmParse pp("spinner");

    pp.add("mode", std::string("fixed"));
    pp.add("turbine", std::string("WTG01"));
    pp.add("hub_debug", false);
    pp.add("inner_prism_theta0", 90.0);
    pp.add("inner_prism_rotrate", 3.5);
    pp.add("inner_prism_azimuth", 15.2);
    pp.add("outer_prism_theta0", 90.0);
    pp.add("outer_prism_rotrate", 6.5);
    pp.add("outer_prism_azimuth", 15.2);
    pp.addarr("lidar_center", amrex::Vector<amrex::Real>{630.0, 192.0, 120.0});
    pp.add("scan_time", 2.0);
    pp.add("num_samples", 984);
    pp.add("beam_length", 270.0);
    pp.add("beam_points", 432);
    pp.add("fixed_yaw", 0);
    pp.add("fixed_roll", 0);
    pp.add("fixed_tilt", 0);

    amr_wind::sampling::DTUSpinnerSampler spinner(sim());
    spinner.initialize("spinner");
    amr_wind::sampling::DTUSpinnerSampler::SampleLocType locs;
    spinner.sampling_locations(locs);

    ASSERT_EQ(locs.size(), 21600);
}

TEST_F(SamplingTest, radar_sampler)
{
    initialize_mesh();
    amrex::ParmParse pp("radar");

    pp.add("num_points", 512);
    pp.addarr("origin", amrex::Vector<amrex::Real>{1.0, 1.0, 1.0});
    pp.add("sampling_frequency", 85.0);
    pp.add("device_sampling_frequency", 30.0);
    pp.add("radar_cone_angle", 0.25);
    pp.add("radar_quadrature_type", std::string("truncated_normal_halfpower"));
    pp.add("radar_npts_azimuth", 5);
    pp.add("radar_beam_length", 100.0);
    pp.add("angular_speed", 30.0);
    pp.add("sweep_angle", 145.0);
    pp.add("reset_time", 0.0);
    pp.addarr(
        "elevation_angles",
        amrex::Vector<amrex::Real>{0.0, 0.1, 0.2, 0.3, 0.4});
    pp.addarr(
        "axis", amrex::Vector<amrex::Real>{0.707106781, 0.707106781, 0.0});
    pp.addarr("vertical_unit_dir", amrex::Vector<amrex::Real>{0.0, 0.0, 1.0});
    pp.add("debug_print", false);

    amr_wind::sampling::RadarSampler radar(sim());
    radar.initialize("radar");
    amr_wind::sampling::RadarSampler::SampleLocType locs;
    radar.sampling_locations(locs);

    ASSERT_EQ(locs.size(), 193536);
}

TEST_F(SamplingTest, sampling_utils)
{
    namespace vs = amr_wind::vs;
    double toler = 1.0e-10;
    vs::Vector unitx{1.0, 0.0, 0.0};
    vs::Vector unity{0.0, 1.0, 0.0};
    vs::Vector unitz{0.0, 0.0, 1.0};
    vs::Vector nunitx{-1.0, 0.0, 0.0};
    vs::Vector ffn{-0.70710678, -0.70710678, 0.0};
    vs::Vector ffnr{-0.70710678, 0.70710678, 0.0};
    vs::Vector result;
    vs::Vector angles{0.0, 180.0, 0.0};

    result = amr_wind::sampling::sampling_utils::reflect(unity, ffn);
    EXPECT_NEAR(result[0], ffnr[0], toler);
    EXPECT_NEAR(result[1], ffnr[1], toler);
    EXPECT_NEAR(result[2], ffnr[2], toler);

    result = amr_wind::sampling::sampling_utils::rotation(angles, unitx);
    EXPECT_NEAR(result[0], nunitx[0], toler);
    EXPECT_NEAR(result[1], nunitx[1], toler);
    EXPECT_NEAR(result[2], nunitx[2], toler);

    result = amr_wind::sampling::sampling_utils::rotate_euler_vec(
        unity, -90.0, unitx);
    EXPECT_NEAR(result[0], unitz[0], toler);
    EXPECT_NEAR(result[1], unitz[1], toler);
    EXPECT_NEAR(result[2], unitz[2], toler);
}

TEST_F(SamplingTest, quadrature)
{
    double toler = 1e-12;
    namespace vs = amr_wind::vs;
    namespace su = amr_wind::sampling::sampling_utils;
    int ntheta = 5;
    double gammav = 0.25 * M_PI / 180.0;
    std::vector<double> weights = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    vs::Vector sr{0.0, 0.0, 1.0};
    std::vector<vs::Vector> rays = {sr, sr, sr, sr, sr, sr, sr, sr, sr, sr, sr,
                                    sr, sr, sr, sr, sr, sr, sr, sr, sr, sr};

    su::spherical_cap_truncated_normal(
        gammav, ntheta, su::NormalRule::HALFPOWER, rays, weights);

    EXPECT_NEAR(weights[0], 1.1826123083219489e-05, toler);
    EXPECT_NEAR(weights[10], 3.004263660003298e-06, toler);
    EXPECT_NEAR(weights[20], 6.6402168628164281e-07, toler);
}

} // namespace amr_wind_tests
