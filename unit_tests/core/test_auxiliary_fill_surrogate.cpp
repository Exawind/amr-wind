#include "aw_test_utils/MeshTest.H"

namespace amr_wind_tests {

namespace {

void auxiliary_fill_boundary(amr_wind::Field& velocity, const int comp = 0)
{
    const int nlevels = velocity.repo().num_active_levels();
    const int ncomp = velocity.num_comp();

    for (int lev = 0; lev < nlevels; ++lev) {

        for (amrex::OrientationIter oit; oit != nullptr; ++oit) {
            auto ori = oit();
            const amrex::Box& minBox =
                velocity.repo().mesh().boxArray(lev).minimalBox();
            amrex::IntVect plo(minBox.loVect());
            amrex::IntVect phi(minBox.hiVect());
            const int normal = ori.coordDir();
            plo[normal] = ori.isHigh() ? minBox.hiVect()[normal] + 1 : -1;
            phi[normal] = ori.isHigh() ? minBox.hiVect()[normal] + 1 : -1;
            const amrex::Box domain_bdy_bx(plo, phi);

            for (amrex::MFIter mfi(velocity(lev), amrex::TilingIfNotGPU());
                 mfi.isValid(); ++mfi) {

                auto sbx = mfi.growntilebox(1);
                const auto& field_location_vector = sbx.type();
                if (!sbx.cellCentered()) {
                    sbx.enclosedCells();
                }
                auto bx = sbx & domain_bdy_bx;
                if (bx.isEmpty()) {
                    continue;
                }

                if (ori.isHigh() &&
                    field_location_vector[ori.coordDir()] == 1) {
                    bx.shift(field_location_vector);
                }

                const auto& dest = velocity(lev).array(mfi);
                amrex::ParallelFor(
                    bx, ncomp,
                    [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
                        dest(i, j, k, n) = static_cast<amrex::Real>(comp + n + 1);
                    });
            }
        }
    }
}

amrex::Real get_field_err(
    amr_wind::Field& field, const bool check_all_ghost, const int comp = 0)
{
    const int lev = 0;
    amrex::Real error_total = 0;
    const int ncomp = field.num_comp();

    error_total += amrex::ReduceSum(
        field(lev), field(lev).nGrow(),
        [=] AMREX_GPU_HOST_DEVICE(
            amrex::Box const& bx,
            amrex::Array4<amrex::Real const> const& f_arr) -> amrex::Real {
            amrex::Real error = 0;

            amrex::Loop(bx, [=, &error](int i, int j, int k) noexcept {
                // Do indices manually to check against box operations in code
                if (ncomp == 1) {
                    // Checking a MAC velocity
                    if (comp == 0) {
                        // Checking u
                        if (check_all_ghost) {
                            // Ghost cells in boundaries
                            if (i == -1 || j == -1 || k == -1 || i == 9 ||
                                j == 8 || k == 8) {
                                error += std::abs(f_arr(i, j, k) - 1.0);
                            } else {
                                error += std::abs(f_arr(i, j, k) - 0.0);
                            }
                        } else {
                            // Valid face-normal cells
                            if ((i == 0 || i == 8) && (j >= 0 && j <= 7) &&
                                (k >= 0 && k <= 7)) {
                                error += std::abs(f_arr(i, j, k) - 1.0);
                            } else {
                                error += std::abs(f_arr(i, j, k) - 0.0);
                            }
                        }
                    } else if (comp == 1) {
                        // Checking v
                        if (check_all_ghost) {
                            // Ghost cells in boundaries
                            if (i == -1 || j == -1 || k == -1 || i == 8 ||
                                j == 9 || k == 8) {
                                error += std::abs(f_arr(i, j, k) - 2.0);
                            } else {
                                error += std::abs(f_arr(i, j, k) - 0.0);
                            }
                        } else {
                            // Valid face-normal cells
                            if ((j == 0 || j == 8) && (i >= 0 && i <= 7) &&
                                (k >= 0 && k <= 7)) {
                                error += std::abs(f_arr(i, j, k) - 2.0);
                            } else {
                                error += std::abs(f_arr(i, j, k) - 0.0);
                            }
                        }
                    } else {
                        // Checking w
                        if (check_all_ghost) {
                            // Ghost cells in boundaries
                            if (i == -1 || j == -1 || k == -1 || i == 8 ||
                                j == 8 || k == 9) {
                                error += std::abs(f_arr(i, j, k) - 3.0);
                            } else {
                                error += std::abs(f_arr(i, j, k) - 0.0);
                            }
                        } else {
                            // Valid face-normal cells
                            if ((k == 0 || k == 8) && (i >= 0 && i <= 7) &&
                                (j >= 0 && j <= 7)) {
                                error += std::abs(f_arr(i, j, k) - 3.0);
                            } else {
                                error += std::abs(f_arr(i, j, k) - 0.0);
                            }
                        }
                    }
                } else {
                    // Checking cell-centered velocity
                    if (check_all_ghost) {
                        // Ghost cells in boundaries
                        if (i <= -1 || j <= -1 || k <= -1 || i >= 8 || j >= 8 ||
                            k >= 8) {
                            error += std::abs(f_arr(i, j, k, 0) - 1.0);
                            error += std::abs(f_arr(i, j, k, 1) - 2.0);
                            error += std::abs(f_arr(i, j, k, 2) - 3.0);
                        } else {
                            error += std::abs(f_arr(i, j, k, 0) - 0.0);
                            error += std::abs(f_arr(i, j, k, 1) - 0.0);
                            error += std::abs(f_arr(i, j, k, 2) - 0.0);
                        }
                    } else {
                        // First ghost cells only
                        const bool is_ibdy =
                            ((i == -1 || i == 8) && (j >= 0 && j <= 7) &&
                             (k >= 0 && k <= 7));
                        const bool is_jbdy =
                            ((j == -1 || j == 8) && (i >= 0 && i <= 7) &&
                             (k >= 0 && k <= 7));
                        const bool is_kbdy =
                            ((k == -1 || k == 8) && (i >= 0 && i <= 7) &&
                             (j >= 0 && j <= 7));
                        if (is_ibdy || is_jbdy || is_kbdy) {
                            error += std::abs(f_arr(i, j, k, 0) - 1.0);
                            error += std::abs(f_arr(i, j, k, 1) - 2.0);
                            error += std::abs(f_arr(i, j, k, 2) - 3.0);
                        } else {
                            error += std::abs(f_arr(i, j, k, 0) - 0.0);
                            error += std::abs(f_arr(i, j, k, 1) - 0.0);
                            error += std::abs(f_arr(i, j, k, 2) - 0.0);
                        }
                    }
                }
            });

            return error;
        });
    amrex::ParallelDescriptor::ReduceRealSum(error_total);
    return error_total;
}
} // namespace

class AuxiliaryFillTest : public MeshTest
{
public:
    void set_up_fields()
    {
        auto& frepo = mesh().field_repo();
        m_vel = &frepo.declare_field("velocity", 3, 2, 1);
        auto mac_vels = frepo.declare_face_normal_field(
            {"u_mac", "v_mac", "w_mac"}, 1, 1, 1);
        m_umac = mac_vels[0];
        m_vmac = mac_vels[1];
        m_wmac = mac_vels[2];

        (*m_vel).setVal(0.);
        (*m_umac).setVal(0.);
        (*m_vmac).setVal(0.);
        (*m_wmac).setVal(0.);
    }

    void prep_test()
    {
        // Default dimensions are n_cell = 8 x 8 x 8
        MeshTest::populate_parameters();
        {
            amrex::Vector<int> periodic{{0, 0, 0}};
            amrex::ParmParse pp("geometry");
            pp.addarr("is_periodic", periodic);
        }
        initialize_mesh();
        set_up_fields();
    }

    amr_wind::Field* m_vel;
    amr_wind::Field* m_umac;
    amr_wind::Field* m_vmac;
    amr_wind::Field* m_wmac;
};

TEST_F(AuxiliaryFillTest, velocity_cc)
{
    prep_test();

    // Test fillpatch and check ghost cells
    auxiliary_fill_boundary(*m_vel);
    const auto err = get_field_err(*m_vel, false);
    EXPECT_DOUBLE_EQ(err, 0.);
}

TEST_F(AuxiliaryFillTest, velocity_xface)
{
    prep_test();

    // Test fillpatch and check ghost cells
    auxiliary_fill_boundary(*m_umac, 0);
    const auto err = get_field_err(*m_umac, true, 0);
    EXPECT_DOUBLE_EQ(err, 0.);
}

TEST_F(AuxiliaryFillTest, velocity_yface)
{
    prep_test();

    // Test fillpatch and check ghost cells
    auxiliary_fill_boundary(*m_vmac, 1);
    const auto err = get_field_err(*m_vmac, true, 1);
    EXPECT_DOUBLE_EQ(err, 0.);
}

TEST_F(AuxiliaryFillTest, velocity_zface)
{
    prep_test();

    // Test fillpatch and check ghost cells
    auxiliary_fill_boundary(*m_wmac, 2);
    const auto err = get_field_err(*m_wmac, true, 2);
    EXPECT_DOUBLE_EQ(err, 0.);
}

} // namespace amr_wind_tests
