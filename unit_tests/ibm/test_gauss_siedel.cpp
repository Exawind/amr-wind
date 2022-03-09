#include "aw_test_utils/AmrexTest.H"
#include "AMReX_ParmParse.H"
#include "AMReX_Gpu.H"
#include "../../amr-wind/immersed_boundary/bluff_body/ghost_cell.H"
#include <AMReX_Array.H>
#include <cmath>
#include <algorithm>
#include <AMReX_Vector.H>

namespace amr_wind_tests {

namespace {

// test matrix A
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
    const amrex::Real tol = 1.0e-12;
    const amrex::Real exact_norm =
        3.74165738677; // For the sample 3x3 linear system
    const int max_num_iterations =
        1000; // Max. number of Gauss-Siedel iterations
    amrex::Real x_norm;

    // Initial guess for 'x'
    amrex::Vector<amrex::Real> x{0.0, 0.0, 0.0};

    // 'b' vector
    const amrex::Vector<amrex::Real> b{9.0, 1.0, 0.0};

    amrex::Array2D<amrex::Real, 0, 2, 0, 2> A;

    // Form matrix A
    form_matrix(A);

    // Gauss-Siedel method to compute 'x'
    for (int i = 0; i < max_num_iterations; i++) {
        amr_wind::ib::gauss_siedel_method(A, b, x);

        x_norm = L2_Norm(x);

        if (amrex::Math::abs(x_norm - exact_norm) < tol) {
            break;
        }
    }

    // Verify if the matrix inversion is correct
    EXPECT_NEAR(x_norm, exact_norm, tol);
}

} // namespace amr_wind_tests
