#include "aw_test_utils/MeshTest.H"
#include "aw_test_utils/iter_tools.H"
#include "aw_test_utils/test_utils.H"
#include "amr-wind/ocean_waves/relaxation_zones/wave_utils_K.H"

namespace amr_wind_tests {

class OceanWavesOpTest : public MeshTest
{
protected:
    void populate_parameters() override
    {
        MeshTest::populate_parameters();

        {
            amrex::ParmParse pp("amr");
            amrex::Vector<int> ncell{{32, 4, 4}};
            pp.add("max_level", 0);
            pp.add("max_grid_size", 4);
            pp.addarr("n_cell", ncell);
        }
        {
            amrex::ParmParse pp("geometry");
            amrex::Vector<amrex::Real> problo{{0.0, 0.0, 0.0}};
            amrex::Vector<amrex::Real> probhi{{10.0, 1.0, 1.0}};

            pp.addarr("prob_lo", problo);
            pp.addarr("prob_hi", probhi);
        }
    }
};

namespace {

void initialize_relaxation_zone_field(
    const amrex::Box& bx, const amrex::Array4<amrex::Real>& targ_farr)
{
    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        targ_farr(i, j, k) = 1.0;
    });
}

void init_relaxation_field(amr_wind::Field& target_field)
{

    run_algorithm(target_field, [&](const int lev, const amrex::MFIter& mfi) {
        auto target_field_arr = target_field(lev).array(mfi);
        const auto& bx = mfi.validbox();
        initialize_relaxation_zone_field(bx, target_field_arr);
    });
}

amrex::Real relaxation_zone_error(amr_wind::Field& comp)
{
    amrex::Real error_total = 0.0;

    for (int lev = 0; lev < comp.repo().num_active_levels(); ++lev) {
    }

    return error_total;
}

} // namespace

TEST_F(OceanWavesOpTest, relaxation_zone)
{
    constexpr double tol = 1.0e-12;

    populate_parameters();
    {
        amrex::ParmParse pp("geometry");
        amrex::Vector<int> periodic{{0, 0, 0}};
        pp.addarr("is_periodic", periodic);
    }

    initialize_mesh();

    auto& repo = sim().repo();
    const int ncomp = 1;
    const int nghost = 3;
    auto& comp_field = repo.declare_field("comp_field", ncomp, nghost);
    auto& target_field = repo.declare_field("target_field", ncomp, nghost);
    comp_field.setVal(0.0);
    target_field.setVal(1.0);
    init_relaxation_field(target_field);
    amrex::Real error_total = relaxation_zone_error(comp_field);
    EXPECT_NEAR(error_total, 0.0, tol);
}

} // namespace amr_wind_tests
