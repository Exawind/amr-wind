#ifndef GHOST_CELL_H_
#define GHOST_CELL_H_

#include <AMReX_FArrayBox.H>
#include <cmath>

namespace amr_wind {
namespace ib {

/** VANDERMONDE MATRIX
 *
 *  phi(x1,x2,x3) = C1*x1*x2*x3 + C2*x1*x2 + C3*x2*x3 + C4*x1*x3 + C5*x1 + C6*x2
 * + C7*x3 + C8
 *
 *  THE Vandemonde Matrix takes the form (underscores denote the cell-centre
 * points)
 *
 *        [ x1*x2*x3|_1     x1*x2|_1    x2*x3|_1 .....  x3|_1      1  ]
 *        | x1*x2*x3|_1     x1*x2|_1    x2*x3|_1 .....  x3|_1      1  |
 *  [V] = | x1*x2*x3|_1     x1*x2|_1    x2*x3|_1 .....  x3|_1      1  |
 *        |   ...             ...         ...            ...      ... |
 *        [ x1*x2*x3|_1     x1*x2|_1    x2*x3|_1 .....  x3|_1      1  ]
 *
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void vandermonde_matrix(
    amrex::Real xIP,
    amrex::Real yIP,
    amrex::Real zIP,
    amrex::Real dx,
    amrex::Real dy,
    amrex::Real dz,
    amrex::Array2D<amrex::Real, 0, 7, 0, 7>& Vmatrix)
{
    int row = 0;

    for (int i = -1; i <= 1; i++) {
        for (int j = -1; j <= 1; j++) {
            for (int k = -1; k <= 1; k++) {
                if (i == 0 || j == 0 || k == 0) continue;

                amrex::Real x1 = xIP + 0.5 * i * dx;
                amrex::Real x2 = yIP + 0.5 * j * dy;
                amrex::Real x3 = zIP + 0.5 * k * dz;

                Vmatrix(row, 0) = x1 * x2 * x3;
                Vmatrix(row, 1) = x1 * x2;
                Vmatrix(row, 2) = x1 * x3;
                Vmatrix(row, 3) = x2 * x3;
                Vmatrix(row, 4) = x1;
                Vmatrix(row, 5) = x2;
                Vmatrix(row, 6) = x3;
                Vmatrix(row, 7) = 1.0;

                row++;
            }
        }
    }
}

/* Gauss-Siedel Method
 *
 * The iterative method solves a square system of N=8
 * linear equations and is of the form
 *
 * A x = b
 *
 * Here, A is the coefficient matrix, and x is the unknown.
 *
 */
template <class T, int Nm1>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void gauss_siedel_method(
    const amrex::Array2D<T, 0, Nm1, 0, Nm1>& A,
    const amrex::Vector<amrex::Real>& b,
    amrex::Vector<amrex::Real>& x)
{
    int N = Nm1 + 1; // Number of entries
    amrex::Real sum_c1, sum_c2;

    for (int row = 0; row < N; row++) {
        sum_c1 = 0.0;
        sum_c2 = 0.0;

        amrex::Real inv_Vrr = 1.0 / A(row, row);

        for (int col = 0; col <= (row - 1); col++) {
            sum_c1 = sum_c1 + A(row, col) * x[col];
        }

        for (int col = row + 1; col <= (N - 1); col++) {
            sum_c2 = sum_c2 + A(row, col) * x[col];
        }

        x[row] = inv_Vrr * (b[row] - sum_c1 - sum_c2);
    }
}

} // namespace ib
} // namespace amr_wind

#endif // GHOST_CELL.H
