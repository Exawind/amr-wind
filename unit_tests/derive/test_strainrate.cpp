/** \file test_strainrate.cpp
 *
 *  Unit tests for amr_wind::Strainrate
 */

#include "aw_test_utils/AmrexTest.H"
#include "derive_K.H"
#include "AnalyticalFunctions.H"

namespace amr_wind_tests {

//! Create unique namespace for this test fixture
class StrainrateTest : public AmrexTest
{};

TEST_F(StrainrateTest, interior)
{

    constexpr double tol = 1.0e-12;

    const int n = 5;
    amrex::Box bx{{0, 0, 0}, {n+1, n+1, n+1}};
    QuadraticAnalyticalFunctions func(n,bx);

    amrex::Real idx = 1.0/func.dx_, idy = 1.0/func.dy_, idz = 1.0/func.dz_;

    for(int i = 2; i < n; ++i) {
        for(int j = 2; j < n; ++j) {
            for(int k = 2; k < n; ++k) {
                const amrex::Real sr = amr_wind::strainrate<amr_wind::StencilInterior>(i,j,k,idx,idy,idz,func.vector_.array());
                EXPECT_NEAR(sr, func.strainrate_.array()(i,j,k), tol);
            }
        }
    }

}

TEST_F(StrainrateTest, ilohi)
{

    constexpr double tol = 1.0e-12;

    const int n = 5;
    amrex::Box bx{{0, 0, 0}, {n+1, n+1, n+1}};
    QuadraticAnalyticalFunctions func(n,bx);

    amrex::Real idx = 1.0/func.dx_, idy = 1.0/func.dy_, idz = 1.0/func.dz_;

    int i = 1;
    for(int j = 2; j < n; ++j) {
        for(int k = 2; k < n; ++k) {
            const amrex::Real sr = amr_wind::strainrate<amr_wind::StencilILO>(i,j,k,idx,idy,idz,func.vector_.array());
            EXPECT_NEAR(sr, func.strainrate_.array()(i,j,k), tol);
        }
    }

    i = n;
    for(int j = 2; j < n; ++j) {
        for(int k = 2; k < n; ++k) {
            const amrex::Real sr = amr_wind::strainrate<amr_wind::StencilIHI>(i,j,k,idx,idy,idz,func.vector_.array());
            EXPECT_NEAR(sr, func.strainrate_.array()(i,j,k), tol);
        }
    }

}


TEST_F(StrainrateTest, jlohi)
{

    constexpr double tol = 1.0e-12;

    const int n = 5;
    amrex::Box bx{{0, 0, 0}, {n+1, n+1, n+1}};
    QuadraticAnalyticalFunctions func(n,bx);

    amrex::Real idx = 1.0/func.dx_, idy = 1.0/func.dy_, idz = 1.0/func.dz_;

    int j = 1;
    for(int i = 2; i < n; ++i) {
        for(int k = 2; k < n; ++k) {
            const amrex::Real sr = amr_wind::strainrate<amr_wind::StencilJLO>(i,j,k,idx,idy,idz,func.vector_.array());
            EXPECT_NEAR(sr, func.strainrate_.array()(i,j,k), tol);
        }
    }

    j = n;
    for(int i = 2; i < n; ++i) {
        for(int k = 2; k < n; ++k) {
            const amrex::Real sr = amr_wind::strainrate<amr_wind::StencilJHI>(i,j,k,idx,idy,idz,func.vector_.array());
            EXPECT_NEAR(sr, func.strainrate_.array()(i,j,k), tol);
        }
    }

}

TEST_F(StrainrateTest, klohi)
{

    constexpr double tol = 1.0e-12;

    const int n = 5;
    amrex::Box bx{{0, 0, 0}, {n+1, n+1, n+1}};
    QuadraticAnalyticalFunctions func(n,bx);

    amrex::Real idx = 1.0/func.dx_, idy = 1.0/func.dy_, idz = 1.0/func.dz_;

    int k = 1;
    for(int i = 2; i < n; ++i) {
        for(int j = 2; j < n; ++j) {
            const amrex::Real sr = amr_wind::strainrate<amr_wind::StencilKLO>(i,j,k,idx,idy,idz,func.vector_.array());
            EXPECT_NEAR(sr, func.strainrate_.array()(i,j,k), tol);
        }
    }

    k = n;
    for(int i = 2; i < n; ++i) {
        for(int j = 2; j < n; ++j) {
            const amrex::Real sr = amr_wind::strainrate<amr_wind::StencilKHI>(i,j,k,idx,idy,idz,func.vector_.array());
            EXPECT_NEAR(sr, func.strainrate_.array()(i,j,k), tol);
        }
    }

}

} // namespace amr_wind_tests
