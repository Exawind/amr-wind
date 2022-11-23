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
    const amrex::Box& bx,
    const amrex::Array4<amrex::Real>& theor_farr,
    amrex::Real dx,
    amrex::Real xlo,
    amrex::Real gen_length)
{
    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        const amrex::Real x = xlo + (i + 0.5) * dx;
        amrex::Real xtilde = std::max(std::min(1. - x / gen_length, 1.0), 0.0);
        theor_farr(i, j, k) =
            std::expm1(std::pow(xtilde, 3.5)) / std::expm1(1.0);
    });
}

void init_relaxation_field(amr_wind::Field& theor_field, amrex::Real gen_length)
{
    const auto& geom = theor_field.repo().mesh().Geom();
    run_algorithm(theor_field, [&](const int lev, const amrex::MFIter& mfi) {
        auto theor_field_arr = theor_field(lev).array(mfi);
        const auto& bx = mfi.validbox();
        const auto& dx = geom[lev].CellSizeArray();
        const auto& problo = geom[lev].ProbLoArray();
        initialize_relaxation_zone_field(
            bx, theor_field_arr, dx[0], problo[0], gen_length);
    });
}

void apply_relaxation_zone_field(
    amr_wind::Field& comp, amr_wind::Field& targ, amrex::Real gen_length)
{

    const auto& geom = comp.repo().mesh().Geom();

    for (int lev = 0; lev < comp.repo().num_active_levels(); ++lev) {
        for (amrex::MFIter mfi(comp(lev)); mfi.isValid(); ++mfi) {
            const auto& gbx = mfi.growntilebox(2);
            const auto& dx = geom[lev].CellSizeArray();
            const auto& problo = geom[lev].ProbLoArray();
            const auto& probhi = geom[lev].ProbHiArray();

            auto comp_arr = comp(lev).array(mfi);
            auto targ_arr = targ(lev).array(mfi);

            amrex::ParallelFor(
                gbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    const amrex::Real x = amrex::min(
                        amrex::max(problo[0] + (i + 0.5) * dx[0], problo[0]),
                        probhi[0]);
                    if (x <= problo[0] + gen_length) {
                        const amrex::Real Gamma =
                            amr_wind::ocean_waves::relaxation_zones::
                                Gamma_generate(x - problo[0], gen_length);
                        comp_arr(i, j, k) = targ_arr(i, j, k) * (1. - Gamma) +
                                            comp_arr(i, j, k) * Gamma;
                    }
                });
        }
    }
}

amrex::Real relaxation_zone_error(amr_wind::Field& comp, amr_wind::Field& targ)
{
    amrex::Real error_total = 0.0;

    for (int lev = 0; lev < comp.repo().num_active_levels(); ++lev) {
        error_total += amrex::ReduceSum(
            comp(lev), targ(lev), 0,
            [=] AMREX_GPU_HOST_DEVICE(
                amrex::Box const& bx,
                amrex::Array4<amrex::Real const> const& comp_arr,
                amrex::Array4<amrex::Real const> const& targ_arr)
                -> amrex::Real {
                amrex::Real error = 0.0;

                amrex::Loop(bx, [=, &error](int i, int j, int k) noexcept {
                    error +=
                        amrex::Math::abs(comp_arr(i, j, k) - targ_arr(i, j, k));
                });

                return error;
            });
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
    amrex::Real gen_length = 4.0;
    auto& comp_field = repo.declare_field("comp_field", ncomp, nghost);
    auto& target_field = repo.declare_field("target_field", ncomp, nghost);
    auto& theoretical_field =
        repo.declare_field("theoretical_field", ncomp, nghost);
    comp_field.setVal(0.0);
    target_field.setVal(1.0);
    init_relaxation_field(theoretical_field, gen_length);
    apply_relaxation_zone_field(comp_field, target_field, gen_length);
    amrex::Real error_total =
        relaxation_zone_error(comp_field, theoretical_field);
    EXPECT_NEAR(error_total, 0.0, tol);
}

} // namespace amr_wind_tests
