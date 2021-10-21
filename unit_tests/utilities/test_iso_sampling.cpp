
#include "aw_test_utils/MeshTest.H"

#include "amr-wind/utilities/sampling/IsoSampling.H"
#include "amr-wind/utilities/sampling/SamplingContainer.H"
#include "amr-wind/utilities/sampling/IsoLineSampler.H"

namespace amr_wind_tests {

namespace {

void init_vof(amr_wind::Field& fld, amrex::Real water_level)
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
                const amrex::Real z = problo[2] + (k + offset) * dx[2];

                for (int d = 0; d < ncomp; d++) {
                    if (z > water_level) {
                        farr(i, j, k, d) = 1.0;
                    } else {
                        farr(i, j, k, d) = 0.0;
                    }
                }
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
};

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
            amrex::Vector<amrex::Real> problo{{0.0, 0.0, 0.0}};
            amrex::Vector<amrex::Real> probhi{{128.0, 128.0, 128.0}};

            pp.addarr("prob_lo", problo);
            pp.addarr("prob_hi", probhi);
        }
    }
};

TEST_F(IsoSamplingTest, sampling)
{
    initialize_mesh();
    auto& repo = sim().repo();
    auto& vof = repo.declare_field("vof", 1, 2);
    init_vof(vof, 0.5*(128.0 - 0.0));

    {
        amrex::ParmParse pp("isosampling");
        pp.add("output_frequency", 1);
        pp.addarr("labels", amrex::Vector<std::string>{"IL1","IL2"});
    }
    {
        amrex::ParmParse pp("isosampling.IL1");
        pp.add("type", std::string("IsoLineSampler"));
        pp.add("field", std::string("vof"));
        pp.add("num_points", 16);
        pp.addarr("start", amrex::Vector<amrex::Real>{64.0,0.0,0.0});
        pp.addarr("end", amrex::Vector<amrex::Real>{64.0,128.0,0.0});
        pp.addarr("orientation", amrex::Vector<amrex::Real>{0.0,0.0,1.0});
    }
    {
        amrex::ParmParse pp("isosampling.IL2");
        pp.add("type", std::string("IsoLineSampler"));
        pp.add("field", std::string("vof"));
        pp.add("num_points", 16);
        pp.addarr("start", amrex::Vector<amrex::Real>{0.0,0.0,0.0});
        pp.addarr("end", amrex::Vector<amrex::Real>{0.0,128.0,0.0});
        pp.addarr("orientation", amrex::Vector<amrex::Real>{1.0,0.0,1.0});
    }

    // amr_wind::sampling::Sampling probes(sim(), "sampling");
    IsoSamplingImpl probes(sim(), "isosampling");
    probes.initialize();
    // Check variable names
    auto var_names = probes.var_names();
    auto nvar = var_names.size();
    auto pcomp_names = probes.pcomp_names();
    auto ncomp = pcomp_names.size();
    EXPECT_EQ(nvar,2);
    EXPECT_EQ(var_names[0],"vof");
    EXPECT_EQ(var_names[1],"vof");
    EXPECT_EQ(ncomp,16);
    probes.post_advance_work();
}

} // namespace amr_wind_tests
