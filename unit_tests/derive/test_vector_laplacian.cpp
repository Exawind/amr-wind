/*
 *  Unit tests for amr_wind::laplacian
 */

#include "aw_test_utils/AmrexTest.H"
#include "derive_K.H"
#include "AnalyticalFunctions.H"

namespace amr_wind_tests {

class VectorlaplacianTest : public AmrexTest
{};

TEST_F(VectorlaplacianTest, interior)
{
    constexpr double tol = 1.0e-10;
    
    const int n = 5;
    amrex::Box bx{{0, 0, 0}, {n+1, n+1, n+1}};
    QuadraticAnalyticalFunctions func(n,bx);
    
    amrex::FArrayBox vectorlaplace(bx, AMREX_SPACEDIM);   
    vectorlaplace.setVal<amrex::RunOn::Host>(0.0);
    
    amrex::Real idx = 1.0/func.dx_, idy = 1.0/func.dy_, idz = 1.0/func.dz_;
    
    for(int i = 2; i < n; ++i) {
        for(int j = 2; j < n; ++j) {
            for(int k = 2; k < n; ++k) {
                amr_wind::laplacian<amr_wind::StencilInterior>(
                    i,j,k,idx,idy,idz,
                    func.vector_.array(), vectorlaplace.array(), AMREX_SPACEDIM);
                for (int l=0; l < AMREX_SPACEDIM; l++)
                    EXPECT_NEAR(vectorlaplace.array()(i,j,k,l),
                                func.vectorlaplace_.array()(i,j,k,l), tol);
            }
        }
    }
    
}

TEST_F(VectorlaplacianTest, ihilo)
{
    constexpr double tol = 1.0e-10;
    
    const int n = 5;
    amrex::Box bx{{0, 0, 0}, {n+1, n+1, n+1}};
    QuadraticAnalyticalFunctions func(n,bx);
    
    amrex::FArrayBox vectorlaplace(bx, AMREX_SPACEDIM);   
    vectorlaplace.setVal<amrex::RunOn::Host>(0.0);  
    
    amrex::Real idx = 1.0/func.dx_, idy = 1.0/func.dy_, idz = 1.0/func.dz_;

    int i=1;
    for(int j = 2; j < n; ++j) {
      for(int k = 2; k < n; ++k) {
        amr_wind::laplacian<amr_wind::StencilILO>(
            i,j,k,idx,idy,idz,
            func.vector_.array(),vectorlaplace.array(),AMREX_SPACEDIM);
        for (int l=0; l < AMREX_SPACEDIM; l++)
            EXPECT_NEAR(vectorlaplace.array()(i,j,k,l),
                        func.vectorlaplace_.array()(i,j,k,l),
                        tol);
      }
    }

    i=n;
    for(int j = 2; j < n; ++j) {
        for(int k = 2; k < n; ++k) {
          amr_wind::laplacian<amr_wind::StencilIHI>(
              i,j,k,idx,idy,idz,
              func.vector_.array(),vectorlaplace.array(),AMREX_SPACEDIM);
          for (int l=0; l < AMREX_SPACEDIM; l++)
              EXPECT_NEAR(vectorlaplace.array()(i,j,k,l),
                          func.vectorlaplace_.array()(i,j,k,l),
                          tol);
        }
    }
    
}

TEST_F(VectorlaplacianTest, jhilo)
{
    constexpr double tol = 1.0e-10;
    
    const int n = 5;
    amrex::Box bx{{0, 0, 0}, {n+1, n+1, n+1}};
    QuadraticAnalyticalFunctions func(n,bx);
    
    amrex::FArrayBox vectorlaplace(bx, AMREX_SPACEDIM);   
    vectorlaplace.setVal<amrex::RunOn::Host>(0.0);  
    
    amrex::Real idx = 1.0/func.dx_, idy = 1.0/func.dy_, idz = 1.0/func.dz_;

    int j=1;
    for(int i = 2; i < n; ++i) {
        for(int k = 2; k < n; ++k) {
            amr_wind::laplacian<amr_wind::StencilJLO>(
                i,j,k,idx,idy,idz,
                func.vector_.array(),vectorlaplace.array(),AMREX_SPACEDIM);
            for (int l=0; l < AMREX_SPACEDIM; l++)
                EXPECT_NEAR(vectorlaplace.array()(i,j,k,l),
                            func.vectorlaplace_.array()(i,j,k,l),
                            tol);
        }
    }
    
    j=n;
    for(int i = 2; i < n; ++i) {
        for(int k = 2; k < n; ++k) {
          amr_wind::laplacian<amr_wind::StencilJHI>(
              i,j,k,idx,idy,idz,
              func.vector_.array(),vectorlaplace.array(),AMREX_SPACEDIM);
          for (int l=0; l < AMREX_SPACEDIM; l++)
              EXPECT_NEAR(vectorlaplace.array()(i,j,k,l),
                          func.vectorlaplace_.array()(i,j,k,l),
                          tol);
        }
    }
    
}

TEST_F(VectorlaplacianTest, khilo)
{
    constexpr double tol = 1.0e-10;
    
    const int n = 5;
    amrex::Box bx{{0, 0, 0}, {n+1, n+1, n+1}};
    LinearAnalyticalFunctions func(n,bx);
    
    amrex::FArrayBox vectorlaplace(bx, AMREX_SPACEDIM);   
    vectorlaplace.setVal<amrex::RunOn::Host>(0.0);  
    
    amrex::Real idx = 1.0/func.dx_, idy = 1.0/func.dy_, idz = 1.0/func.dz_;

    int k=1;
    for(int i = 2; i < n; ++i) {
      for(int j = 2; j < n; ++j) {
        amr_wind::laplacian<amr_wind::StencilKLO>(
            i,j,k,idx,idy,idz,
            func.vector_.array(),vectorlaplace.array(),AMREX_SPACEDIM);
        for (int l=0; l < AMREX_SPACEDIM; l++)
            EXPECT_NEAR(vectorlaplace.array()(i,j,k,l),
                        func.vectorlaplace_.array()(i,j,k,l),
                        tol);
      }
    }

    k=n;
    for(int i = 2; i < n; ++i) {
        for(int j = 2; j < n; ++j) {
          amr_wind::laplacian<amr_wind::StencilKHI>(
              i,j,k,idx,idy,idz,
              func.vector_.array(),vectorlaplace.array(),AMREX_SPACEDIM);
          for (int l=0; l < AMREX_SPACEDIM; l++)
              EXPECT_NEAR(vectorlaplace.array()(i,j,k,l),
                          func.vectorlaplace_.array()(i,j,k,l),
                          tol);
        }
    }
    
}

} // namespace amr_wind_tests
