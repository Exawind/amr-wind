
#include "aw_test_utils/MeshTest.H"

#include "amr-wind/convection/incflo_godunov_minmod.H"
#include "amr-wind/convection/incflo_godunov_upwind.H"

namespace amr_wind_tests {

namespace {

void init_scalar_increasing(amr_wind::Field& fld, int dir)
{
    const int nlevels = fld.repo().num_active_levels();
    for (int lev = 0; lev < nlevels; ++lev) {

        for (amrex::MFIter mfi(fld(lev)); mfi.isValid(); ++mfi) {
            auto bx = mfi.growntilebox(1);
            const auto& farr = fld(lev).array(mfi);

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                farr(i, j, k, 0) = (amrex::Real)(
                    dir == 0 ? i * i : (dir == 1 ? j * j : k * k));
            });
        }
    }
}

void init_scalar_slopechange(amr_wind::Field& fld, int dir, int center)
{
    const int nlevels = fld.repo().num_active_levels();
    for (int lev = 0; lev < nlevels; ++lev) {

        for (amrex::MFIter mfi(fld(lev)); mfi.isValid(); ++mfi) {
            auto bx = mfi.growntilebox(1);
            const auto& farr = fld(lev).array(mfi);

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                farr(i, j, k, 0) = (amrex::Real)(
                    dir == 0 ? amrex::Math::abs(i - center)
                             : (dir == 1 ? amrex::Math::abs(j - center)
                                         : amrex::Math::abs(k - center)));
            });
        }
    }
}

void init_scalar_uniform(amr_wind::Field& fld, amrex::Real cst)
{
    const int nlevels = fld.repo().num_active_levels();
    for (int lev = 0; lev < nlevels; ++lev) {

        for (amrex::MFIter mfi(fld(lev)); mfi.isValid(); ++mfi) {
            auto bx = mfi.growntilebox(1);
            const auto& farr = fld(lev).array(mfi);

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                farr(i, j, k, 0) = cst;
            });
        }
    }
}

void get_output_upwind(
    amr_wind::Field& fld,
    int ii,
    int jj,
    int kk,
    amrex::Real& Im,
    amrex::Real& Ip)
{
    amrex::Gpu::DeviceVector<amrex::Real> dout(2, 0.0);
    auto* dout_ptr = dout.data();
    const int nlevels = fld.repo().num_active_levels();
    for (int lev = 0; lev < nlevels; ++lev) {

        for (amrex::MFIter mfi(fld(lev)); mfi.isValid(); ++mfi) {
            auto bx = mfi.validbox();
            const auto& farr = fld(lev).const_array(mfi);

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                amrex::Real im_tmp, ip_tmp;
                Godunov_upwind_fpu(i, j, k, 0, im_tmp, ip_tmp, farr);
                if (i == ii && j == jj && k == kk) {
                    dout_ptr[0] = im_tmp;
                    dout_ptr[1] = ip_tmp;
                }
            });
        }
    }
    amrex::Gpu::copy(
        amrex::Gpu::deviceToHost, dout.begin(), dout.begin() + 1, &Im);
    amrex::Gpu::copy(amrex::Gpu::deviceToHost, dout.end() - 1, dout.end(), &Ip);
}

void get_output_minmod(
    amr_wind::Field& fld,
    amr_wind::Field& mac_fld,
    amrex::Real dt,
    int ii,
    int jj,
    int kk,
    int dir,
    amrex::Real& Im,
    amrex::Real& Ip)
{
    const int nlevels = fld.repo().num_active_levels();
    amrex::Gpu::DeviceVector<amrex::Real> dout(2, 0.0);
    auto* dout_ptr = dout.data();
    const auto* pbc = fld.bcrec_device().data();
    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& dx = fld.repo().mesh().Geom(lev).CellSizeArray();
        amrex::Box const& domain = fld.repo().mesh().Geom(lev).Domain();
        const auto dlo = amrex::lbound(domain);
        const auto dhi = amrex::ubound(domain);

        for (amrex::MFIter mfi(fld(lev)); mfi.isValid(); ++mfi) {
            auto bx = mfi.validbox();
            const auto& farr = fld(lev).const_array(mfi);
            const auto& vel_mac = mac_fld(lev).const_array(mfi);

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                amrex::Real im_tmp = 0.0;
                amrex::Real ip_tmp = 0.0;
                switch (dir) {
                case 0: {
                    Godunov_minmod_fpu_x(
                        i, j, k, 0, dt, dx[0], im_tmp, ip_tmp, farr, vel_mac,
                        pbc[0], dlo.x, dhi.x);
                    break;
                }
                case 1: {
                    Godunov_minmod_fpu_y(
                        i, j, k, 0, dt, dx[1], im_tmp, ip_tmp, farr, vel_mac,
                        pbc[0], dlo.y, dhi.y);
                    break;
                }
                case 2: {
                    Godunov_minmod_fpu_z(
                        i, j, k, 0, dt, dx[2], im_tmp, ip_tmp, farr, vel_mac,
                        pbc[0], dlo.z, dhi.z);
                    break;
                }
                }
                if (i == ii && j == jj && k == kk) {
                    dout_ptr[0] = im_tmp;
                    dout_ptr[1] = ip_tmp;
                }
            });
        }
    }
    amrex::Gpu::copy(
        amrex::Gpu::deviceToHost, dout.begin(), dout.begin() + 1, &Im);
    amrex::Gpu::copy(amrex::Gpu::deviceToHost, dout.end() - 1, dout.end(), &Ip);
}

} // namespace

class MFluxSchemeTest : public MeshTest
{
protected:
    void populate_parameters() override
    {
        MeshTest::populate_parameters();

        {
            amrex::ParmParse pp("amr");
            amrex::Vector<int> ncell{{nx, nx, nx}};
            pp.add("max_level", 0);
            pp.add("max_grid_size", nx);
            pp.addarr("n_cell", ncell);
        }
        {
            amrex::ParmParse pp("geometry");
            pp.addarr("prob_lo", problo);
            pp.addarr("prob_hi", probhi);
        }
    }
    // Parameters
    const amrex::Vector<amrex::Real> problo{{0.0, -1.0, 0.0}};
    const amrex::Vector<amrex::Real> probhi{{2.0, +1.0, 2.0}};
    const int nx = 5;
    const amrex::Real tol = 1e-12;
};

TEST_F(MFluxSchemeTest, upwind)
{
    initialize_mesh();
    auto& repo = sim().repo();

    // Initialize variables to be reused
    amrex::Real Im = 0.0;
    amrex::Real Ip = 0.0;

    // Initialize field variable
    auto& sc = repo.declare_field("scalar", 1, 1);

    /* -- Constant field portion of test -- */
    // Set up field
    amrex::Real sc_cst = 1.5;
    init_scalar_uniform(sc, sc_cst);
    // Compute interpolated quantities at each face
    int i = nx / 2;
    int j = nx / 2;
    int k = nx / 2;
    get_output_upwind(sc, i, j, k, Im, Ip);
    // Check values
    EXPECT_NEAR(sc_cst, Im, tol);
    EXPECT_NEAR(sc_cst, Ip, tol);

    /* -- Increasing in slope: x, y, z -- */
    for (int n = 0; n < 3; ++n) {
        // Set up field
        init_scalar_increasing(sc, n);
        // Compute interpolated quantities at each face
        get_output_upwind(sc, i, j, k, Im, Ip);
        // Check values
        EXPECT_NEAR(i * i, Im, tol);
        EXPECT_NEAR(i * i, Ip, tol);
    }

    /* -- Change in sign of slope: x, y, z -- */
    for (int n = 0; n < 3; ++n) {
        // Set up field
        init_scalar_slopechange(sc, n, i);
        // Compute interpolated quantities at each face
        get_output_upwind(sc, i, j, k, Im, Ip);
        // Check values
        EXPECT_NEAR(0, Im, tol);
        EXPECT_NEAR(0, Ip, tol);
    }
}

TEST_F(MFluxSchemeTest, minmod)
{
    initialize_mesh();
    auto& repo = sim().repo();

    // Parameters to be referenced during test
    amrex::Real adv_vel = 1.2;
    amrex::Real dt = 0.25;

    // Initialize variables to be reused
    amrex::Real Im = 0.0;
    amrex::Real Ip = 0.0;

    // Initialize field variable
    auto& sc = repo.declare_field("scalar", 1, 1);

    // Initialize mac velocity
    repo.declare_face_normal_field({"umac", "vmac", "wmac"}, 1, 1, 1);
    auto& umac = repo.get_field("umac");
    auto& vmac = repo.get_field("vmac");
    auto& wmac = repo.get_field("wmac");
    umac.setVal(adv_vel);
    vmac.setVal(adv_vel);
    wmac.setVal(-adv_vel);

    /* -- Constant field portion of test -- */
    // Set up field
    amrex::Real sc_cst = 1.5;
    init_scalar_uniform(sc, sc_cst);
    // Compute interpolated quantities at each face
    int i = nx / 2;
    int j = nx / 2;
    int k = nx / 2;
    get_output_minmod(sc, umac, dt, i, j, k, 0, Im, Ip);
    // Check values
    EXPECT_NEAR(sc_cst, Im, tol);
    EXPECT_NEAR(sc_cst, Ip, tol);

    /* -- Increasing in slope: x, y, z -- */
    // Values for checking
    auto ir = (amrex::Real)i;
    amrex::Real dx = sc.repo().mesh().Geom(0).CellSizeArray()[0];
    amrex::Real slp = (std::pow(ir, 2) - std::pow(ir - 1.0, 2)) / dx;
    amrex::Real val_p = std::pow(ir, 2) - slp * 0.5 * (dt * adv_vel - dx);
    amrex::Real val_n = std::pow(ir, 2) + slp * 0.5 * (dt * adv_vel - dx);
    // Set up field (x)
    init_scalar_increasing(sc, 0);
    // Compute interpolated quantities at each face
    get_output_minmod(sc, umac, dt, i, j, k, 0, Im, Ip);
    EXPECT_NEAR(i * i, Im, tol);
    EXPECT_NEAR(val_p, Ip, tol);
    // Set up field (y)
    init_scalar_increasing(sc, 1);
    // Compute interpolated quantities at each face
    get_output_minmod(sc, vmac, dt, i, j, k, 1, Im, Ip);
    EXPECT_NEAR(i * i, Im, tol);
    EXPECT_NEAR(val_p, Ip, tol);
    // Set up field (z)
    init_scalar_increasing(sc, 2);
    // Compute interpolated quantities at each face
    get_output_minmod(sc, wmac, dt, i, j, k, 2, Im, Ip);
    EXPECT_NEAR(val_n, Im, tol);
    EXPECT_NEAR(i * i, Ip, tol);

    /* -- Change in sign of slope -- */
    // Set up field (x)
    init_scalar_slopechange(sc, 0, i);
    // Compute interpolated quantities at each face
    get_output_minmod(sc, umac, dt, i, j, k, 0, Im, Ip);
    // Check values
    EXPECT_NEAR(0, Im, tol);
    EXPECT_NEAR(0, Ip, tol);
}

TEST_F(MFluxSchemeTest, minmodbdy)
{
    initialize_mesh();
    auto& repo = sim().repo();

    // Parameters to be referenced during test
    amrex::Real adv_vel = 1.2;
    amrex::Real dt = 0.25;

    // Initialize variables to be reused
    amrex::Real Im = 0.0;
    amrex::Real Ip = 0.0;

    // Initialize field variable
    auto& sc = repo.declare_field("scalar", 1, 1);

    // Initialize mac velocity
    repo.declare_face_normal_field({"umac", "vmac", "wmac"}, 1, 1, 1);
    auto& umac = repo.get_field("umac");
    auto& vmac = repo.get_field("vmac");
    umac.setVal(-adv_vel);
    vmac.setVal(adv_vel);

    // Initialize BCs
    sc.set_default_fillpatch_bc(sim().time());
    // Set BCs of interest to hoextrap
    sc.bcrec()[0].setLo(0, amrex::BCType::hoextrap);
    sc.bcrec()[0].setHi(1, amrex::BCType::hoextrap);
    // Copy to device to be used
    sc.copy_bc_to_device();

    {
        // Look at lo boundary behavior
        int i = 0;
        int j = nx / 2;
        int k = nx / 2;
        // Values for checking
        auto ir = (amrex::Real)i;
        amrex::Real dx = sc.repo().mesh().Geom(0).CellSizeArray()[0];
        amrex::Real slp = (std::pow(ir, 2) - std::pow(ir - 1.0, 2)) / dx;
        amrex::Real val_n = std::pow(ir, 2) + slp * 0.5 * (dt * adv_vel - dx);
        // Set up field
        init_scalar_increasing(sc, 0);
        // Compute interpolated quantities at each face
        get_output_minmod(sc, umac, dt, i, j, k, 0, Im, Ip);
        EXPECT_NEAR(val_n, Im, tol);
        EXPECT_NEAR(ir * ir, Ip, tol);
    }
    {
        // Look at hi boundary behavior
        int i = nx / 2;
        int j = nx - 1;
        int k = nx / 2;
        // Values for checking
        auto ir = (amrex::Real)j;
        amrex::Real dx = sc.repo().mesh().Geom(0).CellSizeArray()[0];
        amrex::Real slp = (std::pow(ir + 1.0, 2) - std::pow(ir, 2)) / dx;
        amrex::Real val_p = std::pow(ir, 2) - slp * 0.5 * (dt * adv_vel - dx);
        // Set up field
        init_scalar_increasing(sc, 1);
        // Compute interpolated quantities at each face
        get_output_minmod(sc, vmac, dt, i, j, k, 1, Im, Ip);
        EXPECT_NEAR(ir * ir, Im, tol);
        EXPECT_NEAR(val_p, Ip, tol);
    }
}

} // namespace amr_wind_tests
