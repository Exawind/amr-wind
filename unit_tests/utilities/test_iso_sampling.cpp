
#include "aw_test_utils/MeshTest.H"

#include "amr-wind/utilities/sampling/IsoSampling.H"
#include "amr-wind/utilities/sampling/SamplingContainer.H"
#include "amr-wind/utilities/sampling/IsoLineSampler.H"

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
    amrex::Real liwl = -1;

    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& dx = mesh.Geom(lev).CellSizeArray();
        const auto& problo = mesh.Geom(lev).ProbLoArray();

        for (amrex::MFIter mfi(vof_fld(lev)); mfi.isValid(); ++mfi) {
            auto bx = mfi.growntilebox();
            const auto& farr = vof_fld(lev).array(mfi);

            amrex::ParallelFor(bx, [=, &liwl] AMREX_GPU_DEVICE(int i, int j, int k) {
                const amrex::Real z = problo[2] + (k + offset) * dx[2];

                amrex::Real local_vof = std::min(1.0, std::max(0.0,
                    (water_level - (z - offset*dx[2])) / dx[2]));
                farr(i, j, k, d) = local_vof;

                // Only needs to be calculated once
                if (liwl < 0) {
                    // Check for multiphase cell
                    if (local_vof > 0 && local_vof < 1) {
                        liwl = z + (0.5 - local_vof)/(0.0 - local_vof)*dx[2];
                    }
                }
            });
        }
    }
    // If only single-phase cells, no interpolation needed
    if (liwl < 0) liwl = water_level;
    // Return result
    return liwl;
}

class IsoSamplingImpl : public amr_wind::sampling::IsoSampling
{
public:
    IsoSamplingImpl(amr_wind::CFDSim& sim, const std::string& label)
        : amr_wind::sampling::IsoSampling(sim, label)
    {}
    void check_parr(
        const int& i_begin,
        const int& i_end,
        const int& sid,
        const std::string& op,
        amrex::Real* carr,
        const amrex::AmrCore& mesh);
    void check_pos(
        const int& i_begin,
        const int& i_end,
        const int& sid,
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

void IsoSamplingImpl::check_parr(const int& i_begin, const int& i_end,
    const int& sid, const std::string& op, amrex::Real* carr,
    const amrex::AmrCore& mesh)
{
    auto* scont = &(this->sampling_container());
    const int nlevels = mesh.finestLevel() + 1;

    for (int lev = 0; lev < nlevels; ++lev) {

        for (amr_wind::sampling::SamplingContainer::ParIterType pti(
                 *scont, lev); pti.isValid(); ++pti) {
            const int np = pti.numParticles();
            auto& pvec = pti.GetArrayOfStructs()();

            for (int i = i_begin; i < i_end+1; ++i) {
                // Get specified real component
                auto* parr = &(pti.GetStructOfArrays().GetRealData(i))[0];
                // Loop through particles
                auto* pstruct = pvec.data();

                amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE(int ip) noexcept {
                    auto& p = pstruct[ip];
                    // Check if current particle is concerned with current field
                    if (p.idata(amr_wind::sampling::IIx::sid) != sid) return;
                    // Check against reference with specified operation
                    if (op == "=") {
                        EXPECT_EQ(parr[ip], carr[i - i_begin]);
                    } else {
                        if (op == "<") {
                            EXPECT_LT(parr[ip], carr[i - i_begin]);
                        } else {
                            if (op == "~") {
                                EXPECT_NEAR(parr[ip], carr[i - i_begin], tol);
                            }
                        }
                    }
                });
            }
        }
    }
}

void IsoSamplingImpl::check_pos(const int& i_begin, const int& i_end,
    const int& sid, amrex::Real* carr, const amrex::AmrCore& mesh)
{
    auto* scont = &(this->sampling_container());
    const int nlevels = mesh.finestLevel() + 1;

    for (int lev = 0; lev < nlevels; ++lev) {

        for (amr_wind::sampling::SamplingContainer::ParIterType pti(
                 *scont, lev); pti.isValid(); ++pti) {
            const int np = pti.numParticles();
            auto& pvec = pti.GetArrayOfStructs()();

            for (int i = i_begin; i < i_end+1; ++i) {
                // Loop through particles
                auto* pstruct = pvec.data();

                amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE(int ip) noexcept {
                    auto& p = pstruct[ip];
                    // Check if current particle is concerned with current field
                    if (p.idata(amr_wind::sampling::IIx::sid) != sid) return;
                    // Check position against reference with specified component
                    EXPECT_NEAR(p.pos(i), carr[i - i_begin], tol);
                });
            }
        }
    }
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
        }
    }
    void setup_samplers()
    {
        {
            amrex::ParmParse pp("isosampling");
            pp.add("output_frequency", 1);
            pp.addarr("labels", amrex::Vector<std::string>{"IL1", "IL2"});
        }
        {
            amrex::ParmParse pp("isosampling.IL1");
            pp.add("type", std::string("IsoLineSampler"));
            pp.add("field", std::string("vof"));
            pp.add("field_value", 0.5);
            pp.add("num_points", 3);
            pp.addarr("start", IL1_start);
            pp.addarr(
                "end", amrex::Vector<amrex::Real>{
                           probhi[0], IL1_start[1], IL1_start[2]});
            pp.addarr("orientation", amrex::Vector<amrex::Real>{0.0, 0.0, 1.0});
        }
        {
            amrex::ParmParse pp("isosampling.IL2");
            pp.add("type", std::string("IsoLineSampler"));
            pp.add("field", std::string("vof"));
            pp.add("field_value", 0.5);
            pp.add("num_points", 3);
            pp.addarr("start", IL2_start);
            pp.addarr(
                "end", amrex::Vector<amrex::Real>{
                           probhi[0], IL2_start[1], IL2_start[2]});
            pp.addarr("orientation", amrex::Vector<amrex::Real>{0.0, 1.0, 1.0});
        }
    }
    // Parameters to reuse
    const amrex::Real water_level0 = 64.0, water_level1 = 31.5;
    const amrex::Vector<amrex::Real> problo{{0.0, 0.0, 0.0}};
    const amrex::Vector<amrex::Real> probhi{{128.0, 128.0, 128.0}};
    const amrex::Vector<amrex::Real> IL1_start{{0.0, 64.0, 0.0}};
    const amrex::Vector<amrex::Real> IL2_start{{0.0, 0.0, 0.0}};
};

TEST_F(IsoSamplingTest, setup)
{
    initialize_mesh();
    auto& repo = sim().repo();
    auto& vof = repo.declare_field("vof", 1, 2);
    setup_samplers();

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
    EXPECT_EQ(ncomp, 16);
    // Check probe location bounds after iso_initbounds
    amrex::Array<amrex::Real,2> check_two;
    amrex::Array<amrex::Real,1> check_one;
    // Sampler 1
    int sid = 0;
    // Left location should match initial location
    check_two[0] = IL1_start[1];
    check_two[1] = IL1_start[2];
    auto* check_ptr = &(check_two)[0];
    probes.check_parr(4+1,4+2,sid,"=",check_ptr,m_sim.mesh());
    // Left value should be vof = 1 (for this case)
    check_one[0] = 1.0;
    check_ptr = &(check_one)[0];
    probes.check_parr(2,2,sid,"=",check_ptr,m_sim.mesh());
    // Right location should be within bounds
    check_one[0] = probhi[2];
    check_ptr = &(check_one)[0];
    probes.check_parr(4+3+2,4+3+2,sid,"<",check_ptr,m_sim.mesh());
    // Right value should be vof = 0 (for this case)
    check_one[0] = 0.0;
    check_ptr = &(check_one)[0];
    probes.check_parr(3,3,sid,"=",check_ptr,m_sim.mesh());
    // Target value should be recorded
    check_one[0] = 0.5;
    check_ptr = &(check_one)[0];
    probes.check_parr(1,1,sid,"=",check_ptr,m_sim.mesh());
    // Sampler 2
    sid = 1;
    // Right location should be within bounds
    check_two[0] = probhi[1];
    check_two[1] = probhi[2];
    check_ptr = &(check_two)[0];
    probes.check_parr(4+3+1,4+3+2,sid,"<",check_ptr,m_sim.mesh());

}

TEST_F(IsoSamplingTest, once)
{
    initialize_mesh();
    auto& repo = sim().repo();
    auto& vof = repo.declare_field("vof", 1, 2);
    setup_samplers();

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
    amrex::Array<amrex::Real,2> check_two;
    amrex::Array<amrex::Real,1> check_one;
    check_one[0] = 0.5;
    auto* check_ptr = &(check_one)[0];
    probes.check_parr(0,0,sid,"~",check_ptr,m_sim.mesh());
    // Current position should be at water_level
    check_one[0] = liwl;
    check_ptr = &(check_one)[0];
    probes.check_pos(2,2,sid,check_ptr,m_sim.mesh());
    // Sampler 2
    sid = 1;
    // Sample value should be near target
    check_one[0] = 0.5;
    check_ptr = &(check_one)[0];
    probes.check_parr(0,0,sid,"~",check_ptr,m_sim.mesh());
    // Current position should be intersect of water_level and orientation
    check_two[0] = liwl;
    check_two[1] = liwl;
    check_ptr = &(check_two)[0];
    probes.check_pos(1,2,sid,check_ptr,m_sim.mesh());
}

TEST_F(IsoSamplingTest, twice)
{
    initialize_mesh();
    auto& repo = sim().repo();
    auto& vof = repo.declare_field("vof", 1, 2);
    setup_samplers();

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
    amrex::Array<amrex::Real,2> check_two;
    amrex::Array<amrex::Real,1> check_one;
    check_one[0] = 0.5;
    auto* check_ptr = &(check_one)[0];
    probes.check_parr(0,0,sid,"~",check_ptr,m_sim.mesh());
    // Current position should be at water_level
    check_one[0] = liwl;
    check_ptr = &(check_one)[0];
    probes.check_pos(2,2,sid,check_ptr,m_sim.mesh());
    // Sampler 2
    sid = 1;
    // Sample value should be near target
    check_one[0] = 0.5;
    check_ptr = &(check_one)[0];
    probes.check_parr(0,0,sid,"~",check_ptr,m_sim.mesh());
    // Current position should be intersect of water_level and orientation
    check_two[0] = liwl;
    check_two[1] = liwl;
    check_ptr = &(check_two)[0];
    probes.check_pos(1,2,sid,check_ptr,m_sim.mesh());

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
    probes.check_parr(0,0,sid,"~",check_ptr,m_sim.mesh());
    // Current position should be at water_level
    check_one[0] = liwl;
    check_ptr = &(check_one)[0];
    probes.check_pos(2,2,sid,check_ptr,m_sim.mesh());
    // Sampler 2
    sid = 1;
    // Sample value should be near target
    check_one[0] = 0.5;
    check_ptr = &(check_one)[0];
    probes.check_parr(0,0,sid,"~",check_ptr,m_sim.mesh());
    // Current position should be intersect of water_level and orientation
    check_two[0] = liwl;
    check_two[1] = liwl;
    check_ptr = &(check_two)[0];
    probes.check_pos(1,2,sid,check_ptr,m_sim.mesh());

}

} // namespace amr_wind_tests
