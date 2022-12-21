
#include "aw_test_utils/MeshTest.H"

#include "amr-wind/utilities/sampling/Sampling.H"
#include "amr-wind/utilities/sampling/SamplingContainer.H"
#include "amr-wind/utilities/sampling/PlaneSampler.H"

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

        for (amrex::MFIter mfi(fld(lev)); mfi.isValid(); ++mfi) {
            auto bx = mfi.growntilebox();
            const auto& farr = fld(lev).array(mfi);

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                const amrex::Real x = problo[0] + (i + offset) * dx[0];
                const amrex::Real y = problo[1] + (j + offset) * dx[1];
                const amrex::Real z = problo[2] + (k + offset) * dx[2];

                for (int d = 0; d < ncomp; d++) {
                    farr(i, j, k, d) = x + y + z;
                }
            });
        }
    }
}

class SamplingImpl : public amr_wind::sampling::Sampling
{
public:
    SamplingImpl(amr_wind::CFDSim& sim, const std::string& label)
        : amr_wind::sampling::Sampling(sim, label)
    {}

protected:
    void prepare_netcdf_file() override {}
    void process_output() override
    {
        // Test buffer populate for GPU runs
        std::vector<double> buf(
            num_total_particles() * var_names().size(), 0.0);
        sampling_container().populate_buffer(buf);
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

#if 0
    // Old host based implementation for reference
    using PType = amr_wind::sampling::SamplingContainer::ParticleType;
    const int dl = (end[2] - begin[2]) / static_cast<amrex::Real>(npts - 1);
    for (int k = 0; k < npts; ++k) {
        PType p;
        const amrex::Real z = (k + 0.5) * dl;

        p.id() = PType::NextID();
        p.cpu() = amrex::ParallelDescriptor::MyProc();
        p.pos(0) = begin[0];
        p.pos(1) = begin[1];
        p.pos(2) = z;

        p.idata(IIx::sid) = 0;
        p.idata(IIx::nid) = k;

        ptile.push_back(p);

        // Initialize the runtime array stuff
        auto& soa = ptile.GetStructOfArrays();
        for (int ii = 0; ii < ncomp; ++ii) {
            soa.GetRealData(ii).push_back(0.0);
        }
    }
    sc.Redistribute();
#endif

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
    init_field(vel);
    init_field(pres);
    init_field(rho);

    {
        amrex::ParmParse pp("sampling");
        pp.add("output_frequency", 1);
        pp.addarr("labels", amrex::Vector<std::string>{"line1"});
        pp.addarr(
            "fields",
            amrex::Vector<std::string>{"density", "pressure", "velocity"});
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
}

TEST_F(SamplingTest, plane_sampler)
{
    initialize_mesh();

    {
        amrex::ParmParse pp("plane");
        pp.addarr("axis1", amrex::Vector<double>{0.0, 1.0, 0.0});
        pp.addarr("axis2", amrex::Vector<double>{0.0, 0.0, 1.0});
        pp.addarr("origin", amrex::Vector<double>{0.0, 0.0, 0.0});
        pp.addarr("num_points", amrex::Vector<int>{3, 3});
        pp.addarr("offsets", amrex::Vector<double>{-1.0, 1.0});
        pp.addarr("normal", amrex::Vector<double>{1.0, 0.0, 0.0});
    }

    amr_wind::sampling::PlaneSampler plane(sim());
    plane.initialize("plane");
    amr_wind::sampling::PlaneSampler::SampleLocType locs;
    plane.sampling_locations(locs);

    ASSERT_EQ(locs.size(), 3 * 3 * 2);
#if 0
    for (amrex::Long i=0; i < locs.size(); ++i) {
        for (int d=0; d < AMREX_SPACEDIM; ++d) {
            std::cerr << locs[i][d] << " ";
        }
        std::cerr << std::endl;
    }
#endif
}

} // namespace amr_wind_tests
