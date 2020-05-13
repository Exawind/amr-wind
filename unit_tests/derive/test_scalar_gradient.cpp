/** 
 *  Unit tests for amr_wind::gradient
 */

#include "aw_test_utils/AmrexTest.H"
#include "derive_K.H"
#include "AnalyticalFunctions.H"

namespace amr_wind_tests {

class ScalarGradientTest : public AmrexTest
{};

TEST_F(ScalarGradientTest, interior)
{
    constexpr double tol = 1.0e-12;
    
    const int n = 5;
    amrex::Box bx{{0, 0, 0}, {n+1, n+1, n+1}};
    QuadraticAnalyticalFunctions func(n,bx);
    
    amrex::FArrayBox scalargradient(bx, AMREX_SPACEDIM);   
    scalargradient.setVal<amrex::RunOn::Host>(0.0);  
    
    amrex::Real idx = 1.0/func.dx_, idy = 1.0/func.dy_, idz = 1.0/func.dz_;
    
    for(int i = 2; i < n; ++i){
        for(int j = 2; j < n; ++j){
            for(int k = 2; k < n; ++k){
                amr_wind::gradient<amr_wind::StencilInterior>(i,j,k,idx,idy,idz,func.scalar_.array(),scalargradient.array(),1);
                EXPECT_NEAR(scalargradient.array()(i,j,k,0), func.scalargrad_.array()(i,j,k,0), tol);
                EXPECT_NEAR(scalargradient.array()(i,j,k,1), func.scalargrad_.array()(i,j,k,1), tol);
                EXPECT_NEAR(scalargradient.array()(i,j,k,2), func.scalargrad_.array()(i,j,k,2), tol);
            }
        }
    }
    
}

TEST_F(ScalarGradientTest, ihilo)
{
    constexpr double tol = 1.0e-12;
    
    const int n = 5;
    amrex::Box bx{{0, 0, 0}, {n+1, n+1, n+1}};
    LinearAnalyticalFunctions func(n,bx);
    
    amrex::FArrayBox scalargradient(bx, AMREX_SPACEDIM);   
    scalargradient.setVal<amrex::RunOn::Host>(0.0);  
    
    amrex::Real idx = 1.0/func.dx_, idy = 1.0/func.dy_, idz = 1.0/func.dz_;

    int i=1;
    for(int j = 1; j <= n; ++j) {
      for(int k = 1; k <= n; ++k) {
        amr_wind::gradient<amr_wind::StencilILO>(i,j,k,idx,idy,idz,func.scalar_.array(),scalargradient.array(),1);
        EXPECT_NEAR(scalargradient.array()(i,j,k,0), func.scalargrad_.array()(i,j,k,0), tol);
      }
    }

    i=n;
    for(int j = 1; j <= n; ++j) {
        for(int k = 1; k <= n; ++k) {
            amr_wind::gradient<amr_wind::StencilIHI>(i,j,k,idx,idy,idz,func.scalar_.array(),scalargradient.array(),1);
            EXPECT_NEAR(scalargradient.array()(i,j,k,0), func.scalargrad_.array()(i,j,k,0), tol);
        }
    }
    
}

TEST_F(ScalarGradientTest, jhilo)
{
    constexpr double tol = 1.0e-12;
    
    const int n = 5;
    amrex::Box bx{{0, 0, 0}, {n+1, n+1, n+1}};
    LinearAnalyticalFunctions func(n,bx);
    
    amrex::FArrayBox scalargradient(bx, AMREX_SPACEDIM);   
    scalargradient.setVal<amrex::RunOn::Host>(0.0);  
    
    amrex::Real idx = 1.0/func.dx_, idy = 1.0/func.dy_, idz = 1.0/func.dz_;

    int j=1;
    for(int i = 1; i <= n; ++i) {
        for(int k = 1; k <= n; ++k) {
            amr_wind::gradient<amr_wind::StencilJLO>(i,j,k,idx,idy,idz,func.scalar_.array(),scalargradient.array(),1);
            EXPECT_NEAR(scalargradient.array()(i,j,k,1), func.scalargrad_.array()(i,j,k,1), tol);
        }
    }
    
    j=n;
    for(int i = 1; i <= n; ++i) {
        for(int k = 1; k <= n; ++k) {
            amr_wind::gradient<amr_wind::StencilJHI>(i,j,k,idx,idy,idz,func.scalar_.array(),scalargradient.array(),1);
            EXPECT_NEAR(scalargradient.array()(i,j,k,1), func.scalargrad_.array()(i,j,k,1), tol);
        }
    }
    
}

TEST_F(ScalarGradientTest, khilo)
{
    constexpr double tol = 1.0e-12;
    
    const int n = 5;
    amrex::Box bx{{0, 0, 0}, {n+1, n+1, n+1}};
    LinearAnalyticalFunctions func(n,bx);
    
    amrex::FArrayBox scalargradient(bx, AMREX_SPACEDIM);   
    scalargradient.setVal<amrex::RunOn::Host>(0.0);  
    
    amrex::Real idx = 1.0/func.dx_, idy = 1.0/func.dy_, idz = 1.0/func.dz_;

    int k=1;
    for(int i = 1; i <= n; ++i) {
      for(int j = 1; j <= n; ++j) {
        amr_wind::gradient<amr_wind::StencilKLO>(i,j,k,idx,idy,idz,func.scalar_.array(),scalargradient.array(),1);
        EXPECT_NEAR(scalargradient.array()(i,j,k,2), func.scalargrad_.array()(i,j,k,2), tol);
      }
    }

    k=n;
    for(int i = 1; i <= n; ++i) {
        for(int j = 1; j <= n; ++j) {
            amr_wind::gradient<amr_wind::StencilKHI>(i,j,k,idx,idy,idz,func.scalar_.array(),scalargradient.array(),1);
            EXPECT_NEAR(scalargradient.array()(i,j,k,2), func.scalargrad_.array()(i,j,k,2), tol);
        }
    }
    
}

} // namespace amr_wind_tests
