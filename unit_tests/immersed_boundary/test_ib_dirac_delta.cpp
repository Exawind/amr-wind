#include "aw_test_utils/MeshTest.H"

#include "amr-wind/immersed_boundary/IBUtils.H"

namespace amr_wind_tests {
namespace {

class IBDiracDeltaTest : public MeshTest
{
protected:
    void populate_parameters() override
    {
        MeshTest::populate_parameters();

        {
            amrex::ParmParse pp("amr");
            amrex::Vector<int> ncell{{32, 32, 32}};
            pp.add("max_level", 0);
            pp.add("max_grid_size", 16);
            pp.addarr("n_cell", ncell);
        }
        {
            amrex::ParmParse pp("geometry");
            amrex::Vector<amrex::Real> problo{{0.0, 0.0, 0.0}};
            amrex::Vector<amrex::Real> probhi{{32.0, 32.0, 32.0}};

            pp.addarr("prob_lo", problo);
            pp.addarr("prob_hi", probhi);
        }
    }

    void fill_lagrangian_markers(
        amrex::Gpu::DeviceVector<amr_wind::vs::Vector> & pos,
        const int num_pos)
    {
        pos.resize(num_pos);

        for (int ip = 0; ip < num_pos; ++ip) {
            amr_wind::vs::Vector tmp{amrex::Random(), amrex::Random(), amrex::Random()};
            for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                tmp[d] = tmp[d] * 32;
            }
            pos[ip] = tmp;
        }
    }
};

} // namespace

TEST_F(IBDiracDeltaTest, ib_utils)
{
    initialize_mesh();
    auto& frepo = mesh().field_repo();
    auto& field = frepo.declare_field("dummy_field", 3, 3, 1);
    const auto& geom = mesh().Geom();
    const int nlevels = mesh().finestLevel() + 1;

    amrex::Gpu::DeviceVector<amr_wind::vs::Vector> lag_markers;
    fill_lagrangian_markers(lag_markers, 11);

    for (int ip = 0; ip < lag_markers.size(); ++ip) {

        amrex::Real sum  = 0.0;
        amrex::Real sum1 = 0.0;
        amrex::Real sum2 = 0.0;
        amrex::Real sum3 = 0.0;

        for (int lev = 0; lev < nlevels; ++lev) {
            const auto& dx = geom[lev].CellSizeArray();
            const auto& problo = geom[lev].ProbLoArray();

            for (amrex::MFIter mfi(field(lev)); mfi.isValid(); ++mfi) {
                const auto& bx = mfi.tilebox();
                const auto& field_arr = field(lev).array(mfi);

                amrex::ParallelFor(bx, [&] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    const amr_wind::vs::Vector ijk_x{
                        problo[0] + (i + 0.5) * dx[0],
                        problo[1] + (j + 0.5) * dx[1],
                        problo[2] + (k + 0.5) * dx[2]};

                    const auto dist_x = ijk_x - lag_markers[ip];

                    const auto dirac_delta_fac =
                        amr_wind::ib::utils::dirac_delta(dist_x, {dx[0], dx[1], dx[2]});

                    sum  += dirac_delta_fac * dx[0]*dx[1]*dx[2];
                    sum1 += dist_x[0] * dirac_delta_fac * dx[0]*dx[1]*dx[2];
                    sum2 += dist_x[1] * dirac_delta_fac * dx[0]*dx[1]*dx[2];
                    sum3 += dist_x[2] * dirac_delta_fac * dx[0]*dx[1]*dx[2];
                });
            }
        }
        amrex::Print() << sum << " " << sum1 << " " << sum2 << " " << sum3 << std::endl;
    }
}

} // namespace amr_wind_tests
