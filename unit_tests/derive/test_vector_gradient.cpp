/*
 *  Unit tests for amr_wind::gradient
 */

#include "aw_test_utils/AmrexTest.H"
#include "derive_K.H"
#include "AnalyticalFunctions.H"

namespace amr_wind_tests {

class VectorGradientTest : public AmrexTest
{};

TEST_F(VectorGradientTest, interior)
{
    constexpr double tol = 1.0e-12;
    
    const int n = 5;
    amrex::Box bx{{0, 0, 0}, {n+1, n+1, n+1}};
    QuadraticAnalyticalFunctions func(n,bx);
    
    amrex::FArrayBox vectorgrad(bx, AMREX_SPACEDIM*AMREX_SPACEDIM);   
    vectorgrad.setVal<amrex::RunOn::Host>(0.0);
    
    amrex::Real idx = 1.0/func.dx_, idy = 1.0/func.dy_, idz = 1.0/func.dz_;
    
    for(int i = 2; i < n; ++i) {
        for(int j = 2; j < n; ++j) {
            for(int k = 2; k < n; ++k) {
                amr_wind::gradient<amr_wind::StencilInterior>(
                    i,j,k,idx,idy,idz,
                    func.vector_.array(), vectorgrad.array(), AMREX_SPACEDIM);
                for (int l=0; l < AMREX_SPACEDIM * AMREX_SPACEDIM; l++)
                    EXPECT_NEAR(vectorgrad.array()(i,j,k,l),
                                func.vectorgrad_.array()(i,j,k,l), tol);
            }
        }
    }
    
}

TEST_F(VectorGradientTest, ihilo)
{
    constexpr double tol = 1.0e-12;
    
    const int n = 5;
    amrex::Box bx{{0, 0, 0}, {n+1, n+1, n+1}};
    LinearAnalyticalFunctions func(n,bx);
    
    amrex::FArrayBox vectorgrad(bx, AMREX_SPACEDIM*AMREX_SPACEDIM);   
    vectorgrad.setVal<amrex::RunOn::Host>(0.0);  
    
    amrex::Real idx = 1.0/func.dx_, idy = 1.0/func.dy_, idz = 1.0/func.dz_;

    int i=1;
    for(int j = 1; j <= n; ++j) {
      for(int k = 1; k <= n; ++k) {
        amr_wind::gradient<amr_wind::StencilILO>(
            i,j,k,idx,idy,idz,
            func.vector_.array(),vectorgrad.array(),AMREX_SPACEDIM);
        for (int l=0; l < AMREX_SPACEDIM; l++)
            EXPECT_NEAR(vectorgrad.array()(i,j,k,l*AMREX_SPACEDIM+0),
                        func.vectorgrad_.array()(i,j,k,l*AMREX_SPACEDIM+0),
                        tol);
      }
    }

    i=n;
    for(int j = 1; j <= n; ++j) {
        for(int k = 1; k <= n; ++k) {
          amr_wind::gradient<amr_wind::StencilIHI>(
              i,j,k,idx,idy,idz,
              func.vector_.array(),vectorgrad.array(),AMREX_SPACEDIM);
          for (int l=0; l < AMREX_SPACEDIM; l++)
              EXPECT_NEAR(vectorgrad.array()(i,j,k,l*AMREX_SPACEDIM+0),
                          func.vectorgrad_.array()(i,j,k,l*AMREX_SPACEDIM+0),
                          tol);
        }
    }
    
}

TEST_F(VectorGradientTest, jhilo)
{
    constexpr double tol = 1.0e-12;
    
    const int n = 5;
    amrex::Box bx{{0, 0, 0}, {n+1, n+1, n+1}};
    LinearAnalyticalFunctions func(n,bx);
    
    amrex::FArrayBox vectorgrad(bx, AMREX_SPACEDIM*AMREX_SPACEDIM);   
    vectorgrad.setVal<amrex::RunOn::Host>(0.0);  
    
    amrex::Real idx = 1.0/func.dx_, idy = 1.0/func.dy_, idz = 1.0/func.dz_;

    int j=1;
    for(int i = 1; i <= n; ++i) {
        for(int k = 1; k <= n; ++k) {
            amr_wind::gradient<amr_wind::StencilJLO>(
                i,j,k,idx,idy,idz,
                func.vector_.array(),vectorgrad.array(),AMREX_SPACEDIM);
            for (int l=0; l < AMREX_SPACEDIM; l++)
                EXPECT_NEAR(vectorgrad.array()(i,j,k,l*AMREX_SPACEDIM+1),
                            func.vectorgrad_.array()(i,j,k,l*AMREX_SPACEDIM+1),
                            tol);
        }
    }
    
    j=n;
    for(int i = 1; i <= n; ++i) {
        for(int k = 1; k <= n; ++k) {
          amr_wind::gradient<amr_wind::StencilJHI>(
              i,j,k,idx,idy,idz,
              func.vector_.array(),vectorgrad.array(),AMREX_SPACEDIM);
          for (int l=0; l < AMREX_SPACEDIM; l++)
              EXPECT_NEAR(vectorgrad.array()(i,j,k,l*AMREX_SPACEDIM+1),
                          func.vectorgrad_.array()(i,j,k,l*AMREX_SPACEDIM+1),
                          tol);
        }
    }
    
}

TEST_F(VectorGradientTest, khilo)
{
    constexpr double tol = 1.0e-12;
    
    const int n = 5;
    amrex::Box bx{{0, 0, 0}, {n+1, n+1, n+1}};
    LinearAnalyticalFunctions func(n,bx);
    
    amrex::FArrayBox vectorgrad(bx, AMREX_SPACEDIM*AMREX_SPACEDIM);   
    vectorgrad.setVal<amrex::RunOn::Host>(0.0);  
    
    amrex::Real idx = 1.0/func.dx_, idy = 1.0/func.dy_, idz = 1.0/func.dz_;

    int k=1;
    for(int i = 1; i <= n; ++i) {
      for(int j = 1; j <= n; ++j) {
        amr_wind::gradient<amr_wind::StencilKLO>(
            i,j,k,idx,idy,idz,
            func.vector_.array(),vectorgrad.array(),AMREX_SPACEDIM);
        for (int l=0; l < AMREX_SPACEDIM; l++)
            EXPECT_NEAR(vectorgrad.array()(i,j,k,l*AMREX_SPACEDIM+2),
                        func.vectorgrad_.array()(i,j,k,l*AMREX_SPACEDIM+2),
                        tol);
      }
    }

    k=n;
    for(int i = 1; i <= n; ++i) {
        for(int j = 1; j <= n; ++j) {
          amr_wind::gradient<amr_wind::StencilKHI>(
              i,j,k,idx,idy,idz,
              func.vector_.array(),vectorgrad.array(),AMREX_SPACEDIM);
          for (int l=0; l < AMREX_SPACEDIM; l++)
              EXPECT_NEAR(vectorgrad.array()(i,j,k,l*AMREX_SPACEDIM+2),
                          func.vectorgrad_.array()(i,j,k,l*AMREX_SPACEDIM+2),
                          tol);
        }
    }
    
}

} // namespace amr_wind_tests
