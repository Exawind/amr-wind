
#include "aw_test_utils/MeshTest.H"

#include "amr-wind/utilities/sampling/IsoSampling.H"

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

void init_field(amr_wind::Field& fld)
{
    const auto& mesh = fld.repo().mesh();
    const int nlevels = fld.repo().num_active_levels();
    const int ncomp = fld.num_comp();

    amrex::Real offset = 0.0;
    if (fld.field_location() == amr_wind::FieldLoc::CELL) offset = 0.5;

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

                for (int d = 0; d < ncomp; d++)
                    farr(i, j, k, d) = z - 0.05 * (x + y);
            });
        }
    }
}

class IsoSamplingImpl : public amr_wind::sampling::IsoSampling
{
public:
    IsoSamplingImpl(amr_wind::CFDSim& sim, const std::string& label)
        : amr_wind::sampling::IsoSampling(sim, label)
    {}
    int check_parr(
        const int& i_begin,
        const int& i_end,
        const int& sid,
        const std::string& op,
        amrex::Real* carr,
        const amrex::AmrCore& mesh);
    int check_parr(
        const int& i_begin,
        const int& i_end,
        const int& sid,
        int* carr,
        const amrex::AmrCore& mesh);
    int check_pos(
        const int& i_begin,
        const int& i_end,
        const int& sid,
        const std::string& op,
        amrex::Real* carr,
        const amrex::AmrCore& mesh);

protected:
    void prepare_netcdf_file() override {}
    void process_output() override
    {
        // Test buffer populate for GPU runs
        std::vector<double> buf(
            num_total_particles() * realcomps_per_particle(), 0.0);
        sampling_container().populate_buffer(buf);
        std::vector<int> ibuf(
            num_total_particles() * intcomps_per_particle(), 0);
        sampling_container().populate_buffer(ibuf);
    }
    const amrex::Real tol = 1e-8;
};

int IsoSamplingImpl::check_parr(
    const int& i_begin,
    const int& i_end,
    const int& sid,
    const std::string& op,
    amrex::Real* carr,
    const amrex::AmrCore& mesh)
{
    auto* scont = &(this->sampling_container());
    const int nlevels = mesh.finestLevel() + 1;

    // Record number of checks though a reduce loop
    int ncheck = amrex::ReduceSum(
        *scont,
        [=] AMREX_GPU_DEVICE(
            const amr_wind::sampling::SamplingContainer::SuperParticleType&
                p) noexcept -> int {
            if (p.idata(amr_wind::sampling::IIx::sid) != sid) return 0;
            return 1;
        });
    amrex::ParallelDescriptor::ReduceIntSum(ncheck);
    // Multiply by number of fields being checked
    int ncheck_tot = ncheck * (1 + i_end - i_begin);

    // Initialize vectors to store particle information
    amrex::Gpu::DeviceVector<amrex::Real> pdvec(ncheck_tot);
    auto* pdarr = pdvec.data();
    amrex::Vector<amrex::Real> pharr(ncheck_tot);

    for (int lev = 0; lev < nlevels; ++lev) {

        for (amr_wind::sampling::SamplingContainer::ParIterType pti(
                 *scont, lev);
             pti.isValid(); ++pti) {
            const int np = pti.numParticles();
            auto& pvec = pti.GetArrayOfStructs()();

            for (int i = i_begin; i < i_end + 1; ++i) {
                // Get specified real component
                auto* parr = &(pti.GetStructOfArrays().GetRealData(i))[0];
                // Loop through particles
                auto* pstruct = pvec.data();

                amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE(int ip) noexcept {
                    auto& p = pstruct[ip];
                    // Check if current particle is concerned with current
                    // field
                    if (p.idata(amr_wind::sampling::IIx::sid) == sid) {
                        // Check should be performed
                        pdarr
                            [(i - i_begin) * ncheck +
                             p.idata(amr_wind::sampling::IIx::nid)] = parr[ip];
                    }
                });
            }
        }
    }
    // Copy from device to host
    amrex::Gpu::copy(
        amrex::Gpu::deviceToHost, pdvec.begin(), pdvec.end(), pharr.begin());
    // Loop through data and perform checks
    int icheck = 0;
    for (int i = i_begin; i < i_end + 1; ++i) {
        for (int n = 0; n < ncheck; ++n) {
            if (op == "=") {
                EXPECT_EQ(pharr[icheck], carr[i - i_begin]);
            } else {
                if (op == "<") {
                    EXPECT_LT(pharr[icheck], carr[i - i_begin]);
                } else {
                    if (op == "~") {
                        EXPECT_NEAR(pharr[icheck], carr[i - i_begin], tol);
                    }
                }
            }
            ++icheck;
        }
    }
    return ncheck_tot;
}

int IsoSamplingImpl::check_parr(
    const int& i_begin,
    const int& i_end,
    const int& sid,
    int* carr,
    const amrex::AmrCore& mesh)
{
    auto* scont = &(this->sampling_container());
    const int nlevels = mesh.finestLevel() + 1;

    // Record number of checks though a reduce loop
    int ncheck = amrex::ReduceSum(
        *scont,
        [=] AMREX_GPU_DEVICE(
            const amr_wind::sampling::SamplingContainer::SuperParticleType&
                p) noexcept -> int {
            if (p.idata(amr_wind::sampling::IIx::sid) != sid) return 0;
            return 1;
        });
    amrex::ParallelDescriptor::ReduceIntSum(ncheck);
    // Multiply by number of fields being checked
    int ncheck_tot = ncheck * (1 + i_end - i_begin);

    // Initialize vectors to store particle information
    amrex::Gpu::DeviceVector<int> pdvec(ncheck_tot);
    auto* pdarr = pdvec.data();
    amrex::Vector<int> pharr(ncheck_tot);

    for (int lev = 0; lev < nlevels; ++lev) {

        for (amr_wind::sampling::SamplingContainer::ParIterType pti(
                 *scont, lev);
             pti.isValid(); ++pti) {
            const int np = pti.numParticles();
            auto& pvec = pti.GetArrayOfStructs()();

            for (int i = i_begin; i < i_end + 1; ++i) {
                // Get specified real component
                auto* parr = &(pti.GetStructOfArrays().GetIntData(i))[0];
                // Loop through particles
                auto* pstruct = pvec.data();

                amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE(int ip) noexcept {
                    auto& p = pstruct[ip];
                    // Check if current particle is concerned with current
                    // field
                    if (p.idata(amr_wind::sampling::IIx::sid) == sid)
                        // Check should be performed
                        pdarr
                            [(i - i_begin) * ncheck +
                             p.idata(amr_wind::sampling::IIx::nid)] = parr[ip];
                });
            }
        }
    }

    // Copy from device to host
    amrex::Gpu::copy(
        amrex::Gpu::deviceToHost, pdvec.begin(), pdvec.end(), pharr.begin());
    // Loop through data and perform checks
    int icheck = 0;
    for (int i = i_begin; i < i_end + 1; ++i) {
        for (int n = 0; n < ncheck; ++n) {
            // Check against reference
            EXPECT_EQ(pharr[icheck], carr[i - i_begin]);
            ++icheck;
        }
    }

    return ncheck_tot;
}

int IsoSamplingImpl::check_pos(
    const int& i_begin,
    const int& i_end,
    const int& sid,
    const std::string& op,
    amrex::Real* carr,
    const amrex::AmrCore& mesh)
{
    auto* scont = &(this->sampling_container());
    const int nlevels = mesh.finestLevel() + 1;

    // Record number of checks though a reduce loop
    int ncheck = amrex::ReduceSum(
        *scont,
        [=] AMREX_GPU_DEVICE(
            const amr_wind::sampling::SamplingContainer::SuperParticleType&
                p) noexcept -> int {
            if (p.idata(amr_wind::sampling::IIx::sid) != sid) return 0;
            return 1;
        });
    amrex::ParallelDescriptor::ReduceIntSum(ncheck);
    // Multiply by number of fields being checked
    int ncheck_tot = ncheck * (1 + i_end - i_begin);

    // Initialize vectors to store particle information
    amrex::Gpu::DeviceVector<amrex::Real> pdvec(ncheck_tot);
    auto* pdarr = pdvec.data();
    amrex::Vector<amrex::Real> pharr(ncheck_tot);

    for (int lev = 0; lev < nlevels; ++lev) {

        for (amr_wind::sampling::SamplingContainer::ParIterType pti(
                 *scont, lev);
             pti.isValid(); ++pti) {
            const int np = pti.numParticles();
            auto& pvec = pti.GetArrayOfStructs()();

            for (int i = i_begin; i < i_end + 1; ++i) {
                // Loop through particles
                auto* pstruct = pvec.data();

                amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE(int ip) noexcept {
                    auto& p = pstruct[ip];
                    // Check if current particle is concerned with current
                    // field
                    if (p.idata(amr_wind::sampling::IIx::sid) == sid) {
                        // Check should be performed
                        pdarr
                            [(i - i_begin) * ncheck +
                             p.idata(amr_wind::sampling::IIx::nid)] = p.pos(i);
                    }
                });
            }
        }
    }
    // Copy from device to host
    amrex::Gpu::copy(
        amrex::Gpu::deviceToHost, pdvec.begin(), pdvec.end(), pharr.begin());
    // Loop through data and perform checks
    int icheck = 0;
    for (int i = i_begin; i < i_end + 1; ++i) {
        for (int n = 0; n < ncheck; ++n) {
            if (op == "=") {
                EXPECT_EQ(pharr[icheck], carr[i - i_begin]);
            } else {
                if (op == "<") {
                    EXPECT_LT(pharr[icheck], carr[i - i_begin]);
                } else {
                    if (op == "~") {
                        EXPECT_NEAR(pharr[icheck], carr[i - i_begin], tol);
                    }
                }
            }
            ++icheck;
        }
    }
    return ncheck_tot;
}

} // namespace

class IsoSamplingTest : public MeshTest
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

            pp.addarr("prob_lo", problo);
            pp.addarr("prob_hi", probhi);
            pp.addarr("is_periodic", amrex::Vector<int>{{1, 1, 0}});
        }
    }
    void setup_line_samplers()
    {
        {
            amrex::ParmParse pp("isosampling");
            pp.add("output_frequency", 1);
            pp.addarr("labels", amrex::Vector<std::string>{"IL1", "IL2"});
        }
        {
            amrex::ParmParse pp("isosampling.IL1");
            pp.add("type", std::string("LineSampler"));
            pp.add("field", std::string("vof"));
            pp.add("field_value", 0.5);
            pp.add("num_points", npts);
            pp.addarr("start", IL1_start);
            pp.addarr(
                "end",
                amrex::Vector<amrex::Real>{
                    probhi[0] *
                        (1.0 - std::numeric_limits<amrex::Real>::epsilon()),
                    IL1_start[1], IL1_start[2]});
            pp.addarr("orientation", amrex::Vector<amrex::Real>{0.0, 0.0, 1.0});
        }
        {
            amrex::ParmParse pp("isosampling.IL2");
            pp.add("type", std::string("LineSampler"));
            pp.add("field", std::string("vof"));
            pp.add("field_value", 0.5);
            pp.add("num_points", npts);
            pp.addarr("start", IL2_start);
            pp.addarr(
                "end",
                amrex::Vector<amrex::Real>{
                    probhi[0] *
                        (1.0 - std::numeric_limits<amrex::Real>::epsilon()),
                    IL2_start[1], IL2_start[2]});
            pp.addarr("orientation", amrex::Vector<amrex::Real>{0.0, 1.0, 1.0});
        }
    }
    void setup_plane_samplers()
    {
        {
            amrex::ParmParse pp("isosampling");
            pp.add("output_frequency", 1);
            pp.addarr("labels", amrex::Vector<std::string>{"IP1", "IP2"});
        }
        {
            amrex::ParmParse pp("isosampling.IP1");
            pp.add("type", std::string("PlaneSampler"));
            pp.add("field", std::string("vof"));
            pp.add("field_value", 0.5);
            pp.addarr("num_points", amrex::Vector<int>{npts, npts});
            pp.addarr("origin", problo);
            pp.addarr(
                "axis1",
                amrex::Vector<amrex::Real>{
                    probhi[0] *
                        (1.0 - std::numeric_limits<amrex::Real>::epsilon()),
                    0.0, 0.0});
            pp.addarr(
                "axis2",
                amrex::Vector<amrex::Real>{
                    0.0,
                    probhi[1] *
                        (1.0 - std::numeric_limits<amrex::Real>::epsilon()),
                    0.0});
            pp.addarr("orientation", amrex::Vector<amrex::Real>{0.0, 0.0, 1.0});
        }
        {
            amrex::ParmParse pp("isosampling.IP2");
            pp.add("type", std::string("PlaneSampler"));
            pp.add("field", std::string("phi"));
            pp.add("field_value", 0.9 * probhi[0]);
            pp.addarr("num_points", amrex::Vector<int>{npts, npts});
            pp.addarr(
                "origin",
                amrex::Vector<amrex::Real>{
                    0.5 * probhi[0], 0.5 * probhi[1], 0.5 * probhi[2]});
            pp.addarr(
                "axis1", amrex::Vector<amrex::Real>{
                             0.2 * probhi[0], 0.0, 0.2 * probhi[2]});
            pp.addarr(
                "axis2", amrex::Vector<amrex::Real>{
                             0.0, 0.2 * probhi[0], 0.2 * probhi[2]});
            pp.addarr(
                "orientation", amrex::Vector<amrex::Real>{-1.0, -1.0, 1.0});
        }
    }
    // Parameters to reuse
    const amrex::Real water_level0 = 64.0, water_level1 = 31.5;
    const amrex::Vector<amrex::Real> problo{{0.0, 0.0, 0.0}};
    const amrex::Vector<amrex::Real> probhi{{128.0, 128.0, 128.0}};
    const amrex::Vector<amrex::Real> IL1_start{{0.0, 64.0, 0.0}};
    const amrex::Vector<amrex::Real> IL2_start{{0.0, 0.0, 0.0}};
    const int npts = 3;
};

TEST_F(IsoSamplingTest, setup)
{
    initialize_mesh();
    auto& repo = sim().repo();
    auto& vof = repo.declare_field("vof", 1, 2);
    setup_line_samplers();

    init_vof(vof, water_level0);
    auto& m_sim = sim();
    IsoSamplingImpl probes(m_sim, "isosampling");
    probes.initialize();

    // Check variable names
    auto var_names = probes.var_names();
    auto nvar = var_names.size();
    auto pcomp_names = probes.pcomp_names();
    auto ncomp = pcomp_names.size();
    EXPECT_EQ(nvar, 2);
    EXPECT_EQ(var_names[0], "vof");
    EXPECT_EQ(var_names[1], "vof");
    EXPECT_EQ(ncomp, 13);
    // Check probe location bounds after iso_initbounds
    amrex::Array<amrex::Real, 2> check_two;
    amrex::Array<amrex::Real, 1> check_one;
    // Sampler 1
    int sid = 0;
    // Left location should match initial location
    check_two[0] = IL1_start[1];
    check_two[1] = IL1_start[2];
    auto* check_ptr = &(check_two)[0];
    int nleftloc =
        probes.check_parr(4 + 1, 4 + 2, sid, "~", check_ptr, m_sim.mesh());
    ASSERT_EQ(nleftloc, npts * 2);
    // Left value should be vof = 1 (for this case)
    check_one[0] = 1.0;
    check_ptr = &(check_one)[0];
    int nleftval = probes.check_parr(2, 2, sid, "=", check_ptr, m_sim.mesh());
    ASSERT_EQ(nleftval, npts);
    // Right location should be within bounds
    check_one[0] = probhi[2];
    check_ptr = &(check_one)[0];
    int nrightloc = probes.check_parr(
        4 + 3 + 2, 4 + 3 + 2, sid, "<", check_ptr, m_sim.mesh());
    ASSERT_EQ(nrightloc, npts);
    // Right value should be vof = 0 (for this case)
    check_one[0] = 0.0;
    check_ptr = &(check_one)[0];
    int nrightval = probes.check_parr(3, 3, sid, "=", check_ptr, m_sim.mesh());
    ASSERT_EQ(nrightval, npts);
    // Target value should be recorded
    check_one[0] = 0.5;
    check_ptr = &(check_one)[0];
    int ntarget = probes.check_parr(1, 1, sid, "=", check_ptr, m_sim.mesh());
    ASSERT_EQ(ntarget, npts);
    // Sampler 2
    sid = 1;
    // Right location should be within bounds
    check_two[0] = probhi[1];
    check_two[1] = probhi[2];
    check_ptr = &(check_two)[0];
    nrightloc = probes.check_parr(
        4 + 3 + 1, 4 + 3 + 2, sid, "<", check_ptr, m_sim.mesh());
    ASSERT_EQ(nrightloc, npts * 2);
}

TEST_F(IsoSamplingTest, once)
{
    initialize_mesh();
    auto& repo = sim().repo();
    auto& vof = repo.declare_field("vof", 1, 2);
    setup_line_samplers();

    amrex::Real liwl = init_vof(vof, water_level0);
    auto& m_sim = sim();
    IsoSamplingImpl probes(m_sim, "isosampling");
    probes.initialize();

    // Perform isosampling
    probes.post_advance_work();

    // Check results (water_level0)
    // Sampler 1
    int sid = 0;
    // Sample value should be near target
    amrex::Array<amrex::Real, 2> check_two;
    amrex::Array<amrex::Real, 1> check_one;
    check_one[0] = 0.5;
    auto* check_ptr = &(check_one)[0];
    int nsample = probes.check_parr(0, 0, sid, "~", check_ptr, m_sim.mesh());
    ASSERT_EQ(nsample, npts);
    // Current position should be at water_level
    check_one[0] = liwl;
    check_ptr = &(check_one)[0];
    int nsamplepos = probes.check_pos(2, 2, sid, "~", check_ptr, m_sim.mesh());
    ASSERT_EQ(nsamplepos, npts);
    // Sampler 2
    sid = 1;
    // Sample value should be near target
    check_one[0] = 0.5;
    check_ptr = &(check_one)[0];
    nsample = probes.check_parr(0, 0, sid, "~", check_ptr, m_sim.mesh());
    ASSERT_EQ(nsample, npts);
    // Current position should be intersect of water_level and orientation
    check_two[0] = liwl;
    check_two[1] = liwl;
    check_ptr = &(check_two)[0];
    nsamplepos = probes.check_pos(1, 2, sid, "~", check_ptr, m_sim.mesh());
    ASSERT_EQ(nsamplepos, npts * 2);
}

TEST_F(IsoSamplingTest, twice)
{
    initialize_mesh();
    auto& repo = sim().repo();
    auto& vof = repo.declare_field("vof", 1, 2);
    setup_line_samplers();

    amrex::Real liwl = init_vof(vof, water_level0);
    auto& m_sim = sim();
    IsoSamplingImpl probes(m_sim, "isosampling");
    probes.initialize();

    // Perform isosampling
    probes.post_advance_work();

    // Check results (water_level0)
    // Sampler 1
    int sid = 0;
    // Sample value should be near target
    amrex::Array<amrex::Real, 2> check_two;
    amrex::Array<amrex::Real, 1> check_one;
    check_one[0] = 0.5;
    auto* check_ptr = &(check_one)[0];
    int nsample = probes.check_parr(0, 0, sid, "~", check_ptr, m_sim.mesh());
    ASSERT_EQ(nsample, npts);
    // Current position should be at water_level
    check_one[0] = liwl;
    check_ptr = &(check_one)[0];
    int nsamplepos = probes.check_pos(2, 2, sid, "~", check_ptr, m_sim.mesh());
    ASSERT_EQ(nsamplepos, npts);
    // Sampler 2
    sid = 1;
    // Sample value should be near target
    check_one[0] = 0.5;
    check_ptr = &(check_one)[0];
    nsample = probes.check_parr(0, 0, sid, "~", check_ptr, m_sim.mesh());
    ASSERT_EQ(nsample, npts);
    // Current position should be intersect of water_level and orientation
    check_two[0] = liwl;
    check_two[1] = liwl;
    check_ptr = &(check_two)[0];
    nsamplepos = probes.check_pos(1, 2, sid, "~", check_ptr, m_sim.mesh());
    ASSERT_EQ(nsamplepos, npts * 2);

    // Change vof distribution
    liwl = init_vof(vof, water_level1);
    // Perform IsoSampling again
    probes.post_advance_work();

    // Check results (water_level1)
    // Sampler 1
    sid = 0;
    // Sample value should be near target
    check_one[0] = 0.5;
    check_ptr = &(check_one)[0];
    nsample = probes.check_parr(0, 0, sid, "~", check_ptr, m_sim.mesh());
    ASSERT_EQ(nsample, npts);
    // Current position should be at water_level
    check_one[0] = liwl;
    check_ptr = &(check_one)[0];
    nsamplepos = probes.check_pos(2, 2, sid, "~", check_ptr, m_sim.mesh());
    ASSERT_EQ(nsamplepos, npts);
    // Sampler 2
    sid = 1;
    // Sample value should be near target
    check_one[0] = 0.5;
    check_ptr = &(check_one)[0];
    nsample = probes.check_parr(0, 0, sid, "~", check_ptr, m_sim.mesh());
    ASSERT_EQ(nsample, npts);
    // Current position should be intersect of water_level and orientation
    check_two[0] = liwl;
    check_two[1] = liwl;
    check_ptr = &(check_two)[0];
    nsamplepos = probes.check_pos(1, 2, sid, "~", check_ptr, m_sim.mesh());
    ASSERT_EQ(nsamplepos, npts * 2);
}

TEST_F(IsoSamplingTest, plane)
{
    initialize_mesh();
    auto& repo = sim().repo();
    auto& vof = repo.declare_field("vof", 1, 2);
    auto& phi = repo.declare_field("phi", 1, 2);
    setup_plane_samplers();

    init_vof(vof, probhi[2] + 1.0); // Put interface above domain
    init_field(phi);
    auto& m_sim = sim();
    IsoSamplingImpl probes(m_sim, "isosampling");
    probes.initialize();

    // Perform isosampling
    probes.post_advance_work();

    // Check results
    // Sampler 1
    int sid = 0;
    // Sampler cannot find target value
    amrex::Array<int, 1> check_int_one;
    check_int_one[0] = -3;
    auto* check_int_ptr = &(check_int_one)[0];
    int niflag = probes.check_parr(0, 0, sid, check_int_ptr, m_sim.mesh());
    ASSERT_EQ(niflag, npts * npts);
    // Sampler 2
    sid = 1;
    // Sampler can find sign change and bisection reaches target
    check_int_one[0] = 2;
    check_int_ptr = &(check_int_one)[0];
    niflag = probes.check_parr(0, 0, sid, check_int_ptr, m_sim.mesh());
    ASSERT_EQ(niflag, npts * npts);
    // Current height should be below highest possible height for target
    amrex::Array<amrex::Real, 1> check_one;
    check_one[0] = 0.9 * (probhi[0] + 0.05 * (probhi[1] + probhi[2]));
    auto* check_ptr = &(check_one)[0];
    int nsamplepos = probes.check_pos(2, 2, sid, "<", check_ptr, m_sim.mesh());
    ASSERT_EQ(nsamplepos, npts * npts);
    // Current value should be near target
    check_one[0] = 0.9 * probhi[0];
    check_ptr = &(check_one)[0];
    int nsample = probes.check_parr(0, 0, sid, "~", check_ptr, m_sim.mesh());
    ASSERT_EQ(nsample, npts * npts);
}

} // namespace amr_wind_tests
