#include "aw_test_utils/MeshTest.H"
#include "amr-wind/fvm/gradient.H"
#include "aw_test_utils/iter_tools.H"

namespace amr_wind_tests {

class FvmOpTest : public MeshTest
{};

namespace {

amrex::Real phi_eval(amrex::Real x, amrex::Real y, amrex::Real z, int index)
{
    amrex::Real phi[3];

    phi[0] = 5.0 * x * y * z + 4.0 * x * x + 3.0 * y * y - 2.3 * z * z +
             1.3 * y * z + 3.8 * x * z + 9.4 * x * y + 3.4 * x + 2.0 * x + 3.14;
    phi[1] = 8.265 * x * y * z + 1.924 * x * x + 0.923 * y * y - 8.65 * z * z +
             2.834 * y * z + 9.812 * x * z + 4.12 * x * y + 1.0921 * x + 131.0;
    phi[2] = 0.2994917334 * x * y * z + 0.8228505242180528 * x * x +
             0.1205690389747357 * y * y - 0.9604194947750178 * z * z +
             0.37039876870122856 * y * z + 0.6241061971479255 * x * z +
             0.7511807179790003 * x * y;

    return phi[index];
}

amrex::Real dphi_eval(amrex::Real x, amrex::Real y, amrex::Real z, int index)
{
    amrex::Real dphi[9];

    dphi[0] = 5.0 * y * z + 8.0 * x + 3.8 * z + 9.4 * y + 3.4 + 2.0;
    dphi[1] = 5.0 * x * z + 6.0 * y + 1.3 * z + 9.4 * x;
    dphi[2] = 5.0 * x * y - 4.6 * z + 1.3 * y + 3.8 * x;
    dphi[3] = 8.265 * y * z + 3.848 * x + 9.812 * z + 4.12 * y + 1.0921;
    dphi[4] = 8.265 * x * z + 1.846 * y + 2.834 * z + 4.12 * x;
    dphi[5] = 8.265 * x * y - 17.3 * z + 2.834 * y + 9.812 * x;

    dphi[6] = 0.2994917334 * y * z + 2.0 * 0.8228505242180528 * x +
              0.6241061971479255 * z + 0.7511807179790003 * y;
    dphi[7] = 0.2994917334 * x * z + 2.0 * 0.1205690389747357 * y +
              0.37039876870122856 * z + 0.7511807179790003 * x;
    dphi[8] = 0.2994917334 * x * y - 2.0 * 0.9604194947750178 * z +
              0.37039876870122856 * y + 0.6241061971479255 * x;

    return dphi[index];
}

void init_vel_field(
    const amrex::Geometry& geom,
    const amrex::Box& bx,
    const amrex::Array4<amrex::Real>& vel)
{
    const auto& problo = geom.ProbLoArray();
    const auto& probhi = geom.ProbHiArray();
    const auto& dx = geom.CellSizeArray();

    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        const amrex::Real x = amrex::min(
            amrex::max(problo[0] + (i + 0.5) * dx[0], problo[0]), probhi[0]);
        const amrex::Real y = amrex::min(
            amrex::max(problo[1] + (j + 0.5) * dx[1], problo[1]), probhi[1]);
        const amrex::Real z = amrex::min(
            amrex::max(problo[2] + (k + 0.5) * dx[2], problo[2]), probhi[2]);

        vel(i, j, k, 0) = phi_eval(x, y, z, 0);
        vel(i, j, k, 1) = phi_eval(x, y, z, 1);
        vel(i, j, k, 2) = phi_eval(x, y, z, 2);
    });
}

} // namespace

double grad_test_impl(amr_wind::Field& vel)
{
    const auto& mesh = vel.repo().mesh();
    const int nghost = vel.num_grow()[0];

    run_algorithm(vel, [&](const int lev, const amrex::MFIter& mfi) {
        const auto bx = amrex::grow(mfi.validbox(), nghost);
        const auto& vel_arr = vel(lev).array(mfi);
        init_vel_field(mesh.Geom(lev), bx, vel_arr);
    });

    auto grad_vel = amr_wind::fvm::gradient(vel);

    amrex::Real error_total = 0.0;

    const auto& problo = mesh.Geom(0).ProbLoArray();
    const auto& dx = mesh.Geom(0).CellSizeArray();

    error_total += amrex::ReduceSum(
        (*grad_vel)(0), vel(0), 0,
        [=] AMREX_GPU_HOST_DEVICE(
            amrex::Box const& bx, amrex::Array4<amrex::Real const> const& gvel,
            amrex::Array4<amrex::Real const> const& vel) -> amrex::Real {
            amrex::Real error = 0.0;

            amrex::Loop(bx, [=, &error](int i, int j, int k) noexcept {
                const amrex::Real x = problo[0] + (i + 0.5) * dx[0];
                const amrex::Real y = problo[1] + (j + 0.5) * dx[1];
                const amrex::Real z = problo[2] + (k + 0.5) * dx[2];

                for (int n = 0; n < 9; ++n) {
                    error += amrex::Math::abs(
                        gvel(i, j, k, n) - dphi_eval(x, y, z, n));
                }
            });

            return error;
        });

    return error_total;
}

TEST_F(FvmOpTest, gradient)
{

    populate_parameters();
    {
        amrex::ParmParse pp("geometry");
        amrex::Vector<int> periodic{{0, 0, 0}};
        pp.addarr("is_periodic", periodic);
    }
    initialize_mesh();

    auto& repo = sim().repo();
    const int nghost = 1;
    auto& vel = repo.declare_field("vel", 3, nghost);

    vel.setVal({666.0, 777.0, 888.0}, nghost);

    auto error_total = grad_test_impl(vel);

    amrex::Print() << "error: " << error_total << std::endl;
}

} // namespace amr_wind_tests
