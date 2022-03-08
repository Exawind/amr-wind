#include "aw_test_utils/AmrexTest.H"
#include "AMReX_ParmParse.H"
#include "AMReX_Gpu.H"
#include "../../amr-wind/immersed_boundary/bluff_body/ghost_cell.H"

namespace amr_wind_tests {

namespace {

void form_matrix(amrex::Array2D<amrex::Real, 0, 2, 0, 2>& A)
{
    A(0, 0) = 1.0;
    A(0, 1) = 1.0;
    A(0, 2) = 2.0;

    A(1, 0) = 2.0;
    A(1, 1) = 4.0;
    A(1, 2) = -3.0;

    A(2, 0) = 3.0;
    A(2, 1) = 6.0;
    A(2, 2) = -5.0;
}

/* Gauss-Siedel Method
 *
 * The iterative method solves a square system of N=3
 * linear equations and is of the form
 *
 * A x = b
 *
 * Here, A is the coefficient matrix, and x is the unknown.
 *
 */
void gauss_siedel_method(
    amrex::Array2D<amrex::Real, 0, 2, 0, 2>& A,
    const amrex::Vector<amrex::Real>& b,
    amrex::Vector<amrex::Real>& x,
    const int NIter)
{
    const int Np = 3;
    int n = 0;

    amrex::Real sum_c1;
    amrex::Real sum_c2;

    while (n < NIter) {
        for (int row = 0; row < Np; row++) {
            sum_c1 = 0.0;
            sum_c2 = 0.0;

            amrex::Real inv_Vrr = 1.0 / A(row, row);

            for (int col = 0; col <= (row - 1); col++) {
                sum_c1 = sum_c1 + A(row, col) * x[col];
            }

            for (int col = row + 1; col <= (Np - 1); col++) {
                sum_c2 = sum_c2 + A(row, col) * x[col];
            }

            x[row] = inv_Vrr * (b[row] - sum_c1 - sum_c2);
        }
        n++;
    }
}

amrex::Real L2_Norm(const amrex::Vector<amrex::Real>& x)
{
    amrex::Real norm = 0.0;

    for (int i = 0; i < x.size(); i++) {
        norm = norm + std::pow(x[i], 2);
    }

    norm = std::sqrt(norm);

    return norm;
}

} // namespace

// Solves 8x8 linear system using Gauss-Siedel method
class GaussSiedelMethodTest : public AmrexTest
{};

TEST_F(GaussSiedelMethodTest, matrix_inverse)
{
    constexpr double tol = 1.0e-3;
    const amrex::Real exact_norm =
        3.74165738677; // For the sample 3x3 linear system

    const int NIter = 100;
    amrex::Real x_norm;

    amrex::Array2D<amrex::Real, 0, 2, 0, 2> A;

    // Initial guess for 'x'
    amrex::Vector<amrex::Real> x{0.0, 0.0, 0.0};

    // 'b' vector
    const amrex::Vector<amrex::Real> b{9.0, 1.0, 0.0};

    // Form matrix A
    form_matrix(A);

    // Gauss-Siedel method to compute 'x'
    gauss_siedel_method(A, b, x, NIter);

    // L2 Norm of 'x'
    x_norm = L2_Norm(x);

    // Verify if the matrix inversion is correct
    EXPECT_NEAR(x_norm, exact_norm, tol);
}

} // namespace amr_wind_tests
