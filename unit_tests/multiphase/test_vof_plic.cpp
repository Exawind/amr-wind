#include "aw_test_utils/MeshTest.H"
#include "aw_test_utils/iter_tools.H"
#include "aw_test_utils/test_utils.H"
#include "amr-wind/equation_systems/vof/volume_fractions.H"

namespace amr_wind_tests {

class VOFOpTest : public MeshTest
{
protected:
    void populate_parameters() override
    {
        MeshTest::populate_parameters();

        {
            amrex::ParmParse pp("amr");
            amrex::Vector<int> ncell{{4, 4, 4}};
            pp.add("max_level", 0);
            pp.add("max_grid_size", 4);
            pp.addarr("n_cell", ncell);
        }
        {
            amrex::ParmParse pp("geometry");
            amrex::Vector<amrex::Real> problo{{0.0, 0.0, 0.0}};
            amrex::Vector<amrex::Real> probhi{{1.0, 1.0, 1.0}};

            pp.addarr("prob_lo", problo);
            pp.addarr("prob_hi", probhi);
        }
    }
};

namespace {

void initialize_volume_fractions(
    const amrex::Geometry&,
    const amrex::Box& bx,
    const int,
    amrex::Array4<amrex::Real>& vof_arr)
{

    // grow the box by 1 so that x,y,z go out of bounds and min(max()) corrects
    // it and it fills the ghosts with wall values
    amrex::ParallelFor(grow(bx, 1), [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        if (i + k > 3) vof_arr(i, j, k) = 0.0;
        if (i + k == 3) vof_arr(i, j, k) = 0.5;
        if (i + k < 3) vof_arr(i, j, k) = 1.0;
    });
}

amrex::Real normal_vector_test_impl(amr_wind::Field& vof, const int pdegree)
{

    auto& geom = vof.repo().mesh().Geom();

    run_algorithm(vof, [&](const int lev, const amrex::MFIter& mfi) {
        auto vof_arr = vof(lev).array(mfi);
        const auto& bx = mfi.validbox();
        initialize_volume_fractions(geom[lev], bx, pdegree, vof_arr);
    });

    amrex::Real error_total = 0.0;

    run_algorithm(vof, [&](const int lev, const amrex::MFIter& mfi) {
        auto vof_arr = vof(lev).array(mfi);
        const auto& bx = mfi.validbox();
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            amrex::Real mx, my, mz;
            amr_wind::multiphase::mixed_youngs_central_normal(
                i, j, k, vof_arr, mx, my, mz);
            amrex::Print() << vof_arr(i, j, k) << " (" << mx << "," << my << ","
                           << mz << ") " << std::endl;
        });
    });
    return error_total;
}

} // namespace

TEST_F(VOFOpTest, interface_normal)
{

    constexpr double tol = 1.0e-11;

    populate_parameters();
    {
        amrex::ParmParse pp("geometry");
        amrex::Vector<int> periodic{{0, 0, 0}};
        pp.addarr("is_periodic", periodic);
    }

    initialize_mesh();

    auto& repo = sim().repo();
    const int ncomp = 3;
    const int nghost = 1;
    auto& vof = repo.declare_field("vof", ncomp, nghost);

    const int pdegree = 2;
    auto error_total = normal_vector_test_impl(vof, pdegree);

    amrex::ParallelDescriptor::ReduceRealSum(error_total);

    EXPECT_NEAR(error_total, 0.0, tol);
}

} // namespace amr_wind_tests
