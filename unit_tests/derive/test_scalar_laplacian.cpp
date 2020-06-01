/** 
 *  Unit tests for amr_wind::laplacian
 */

#include "aw_test_utils/AmrexTest.H"
#include "amr-wind/derive/derive_K.H"
#include "AnalyticalFunctions.H"

namespace amr_wind_tests {

class ScalarlaplacianTest : public AmrexTest
{};

TEST_F(ScalarlaplacianTest, interior)
{
    constexpr double tol = 1.0e-11;
    
    const int n = 5;
    amrex::Box bx{{0, 0, 0}, {n+1, n+1, n+1}};
    QuadraticAnalyticalFunctions func(n,bx);
    
    amrex::FArrayBox scalarlaplacian(bx, 1);   
    scalarlaplacian.setVal<amrex::RunOn::Host>(0.0);  
    
    amrex::Real idx = 1.0/func.dx_, idy = 1.0/func.dy_, idz = 1.0/func.dz_;
    
    for(int i = 2; i < n; ++i){
        for(int j = 2; j < n; ++j){
            for(int k = 2; k < n; ++k){
                amr_wind::laplacian<amr_wind::StencilInterior>(i,j,k,idx,idy,idz,func.scalar_.array(),scalarlaplacian.array(),1);
                EXPECT_NEAR(scalarlaplacian.array()(i,j,k), func.scalarlaplace_.array()(i,j,k), tol);
            }
        }
    }
    
}

TEST_F(ScalarlaplacianTest, ihilo)
{
    constexpr double tol = 1.0e-10;
    
    const int n = 5;
    amrex::Box bx{{0, 0, 0}, {n+1, n+1, n+1}};
    QuadraticAnalyticalFunctions func(n,bx);
    
    amrex::FArrayBox scalarlaplacian(bx, 1);   
    scalarlaplacian.setVal<amrex::RunOn::Host>(0.0);  
    
    amrex::Real idx = 1.0/func.dx_, idy = 1.0/func.dy_, idz = 1.0/func.dz_;

    int i=1;
    for(int j = 2; j < n; ++j) {
      for(int k = 2; k < n; ++k) {
        amr_wind::laplacian<amr_wind::StencilILO>(i,j,k,idx,idy,idz,func.scalar_.array(),scalarlaplacian.array(),1);
        EXPECT_NEAR(scalarlaplacian.array()(i,j,k), func.scalarlaplace_.array()(i,j,k), tol);
      }
    }

    i=n;
    for(int j = 2; j < n; ++j) {
        for(int k = 2; k < n; ++k) {
            amr_wind::laplacian<amr_wind::StencilIHI>(i,j,k,idx,idy,idz,func.scalar_.array(),scalarlaplacian.array(),1);
            EXPECT_NEAR(scalarlaplacian.array()(i,j,k), func.scalarlaplace_.array()(i,j,k), tol);
        }
    }
    
}

TEST_F(ScalarlaplacianTest, jhilo)
{
    constexpr double tol = 1.0e-10;
    
    const int n = 5;
    amrex::Box bx{{0, 0, 0}, {n+1, n+1, n+1}};
    QuadraticAnalyticalFunctions func(n,bx);
    
    amrex::FArrayBox scalarlaplacian(bx, AMREX_SPACEDIM);   
    scalarlaplacian.setVal<amrex::RunOn::Host>(0.0);  
    
    amrex::Real idx = 1.0/func.dx_, idy = 1.0/func.dy_, idz = 1.0/func.dz_;

    int j=1;
    for(int i = 2; i < n; ++i) {
        for(int k = 2; k < n; ++k) {
            amr_wind::laplacian<amr_wind::StencilJLO>(i,j,k,idx,idy,idz,func.scalar_.array(),scalarlaplacian.array(),1);
            EXPECT_NEAR(scalarlaplacian.array()(i,j,k), func.scalarlaplace_.array()(i,j,k), tol);
        }
    }
    
    j=n;
    for(int i = 2; i < n; ++i) {
        for(int k = 2; k < n; ++k) {
            amr_wind::laplacian<amr_wind::StencilJHI>(i,j,k,idx,idy,idz,func.scalar_.array(),scalarlaplacian.array(),1);
            EXPECT_NEAR(scalarlaplacian.array()(i,j,k), func.scalarlaplace_.array()(i,j,k), tol);
        }
    }
    
}

TEST_F(ScalarlaplacianTest, khilo)
{
    constexpr double tol = 1.0e-10;
    
    const int n = 5;
    amrex::Box bx{{0, 0, 0}, {n+1, n+1, n+1}};
    QuadraticAnalyticalFunctions func(n,bx);
    
    amrex::FArrayBox scalarlaplacian(bx, AMREX_SPACEDIM);   
    scalarlaplacian.setVal<amrex::RunOn::Host>(0.0);  
    
    amrex::Real idx = 1.0/func.dx_, idy = 1.0/func.dy_, idz = 1.0/func.dz_;

    int k=1;
    for(int i = 2; i < n; ++i) {
      for(int j = 2; j < n; ++j) {
        amr_wind::laplacian<amr_wind::StencilKLO>(i,j,k,idx,idy,idz,func.scalar_.array(),scalarlaplacian.array(),1);
        EXPECT_NEAR(scalarlaplacian.array()(i,j,k), func.scalarlaplace_.array()(i,j,k), tol);
      }
    }

    k=n;
    for(int i = 2; i < n; ++i) {
        for(int j = 2; j < n; ++j) {
            amr_wind::laplacian<amr_wind::StencilKHI>(i,j,k,idx,idy,idz,func.scalar_.array(),scalarlaplacian.array(),1);
            EXPECT_NEAR(scalarlaplacian.array()(i,j,k), func.scalarlaplace_.array()(i,j,k), tol);
        }
    }
    
}

} // namespace amr_wind_tests
