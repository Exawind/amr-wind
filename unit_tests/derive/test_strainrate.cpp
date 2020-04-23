/** \file test_strainrate.cpp
 *
 *  Unit tests for amr_wind::Strainrate
 */

#include "aw_test_utils/AmrexTest.H"
#include "derive_K.H"

namespace amr_wind_tests {

//! Create unique namespace for this test fixture
class StrainrateTest : public AmrexTest
{};

TEST_F(StrainrateTest, interior)
{

    constexpr double tol = 1.0e-12;
    
    const int n = 4;
    amrex::Box bx{{0, 0, 0}, {n, n, n}};
    amrex::FArrayBox velocity(bx, AMREX_SPACEDIM);
    amrex::FArrayBox strainrate(bx, 1);

    velocity.setVal<amrex::RunOn::Host>(0.0);

    int i = 2, j = 2, k = 2;
    amrex::Real dx = 0.1, dy = 0.2 + 0.01*amrex::Random(), dz = 0.3;
    amrex::Real idx = 1.0/dx, idy = 1.0/dy, idz = 1.0/dz;
    
    
    amrex::Real sr = amr_wind::strainrate<amr_wind::StencilInterior>(i,j,k,idx,idy,idz,velocity.array());
    EXPECT_NEAR(sr, 0.0, tol);
    
    auto vel_arr = velocity.array();
    auto sr_arr = strainrate.array();
    
    // Interior stencil should be exact for quadratic polynomials
    for(i = 0; i <= n; ++i){
        for(j = 0; j <= n; ++j){
            for(k = 0; k <= n; ++k){
                amrex::Real x = 2.0 + i*dx;
                amrex::Real y = 2.5 + j*dy;
                amrex::Real z = 5.2 + k*dz;
                
                vel_arr(i,j,k,0) = 5.0*x*y*z + 4.0*x*x + 3.0*y*y - 2.3*z*z + 1.3*y*z + 3.8*x*z + 9.4*x*y + 3.4*x+ 2.0*x + 3.14;
                amrex::Real ux = 5.0*y*z + 8.0*x + 3.8*z + 9.4*y + 3.4 + 2.0;
                amrex::Real uy = 5.0*x*z + 6.0*y + 1.3*z + 9.4*x;
                amrex::Real uz = 5.0*x*y -4.6*z + 1.3*y + 3.8*x;
                
                vel_arr(i,j,k,1) = 3.0*y*y - 2.3*z*z + 1.3*y*z + 3.8*x*z - 2.9*y + 2.0*x + 2.8;
                amrex::Real vx = 3.8*z + 2.0;
                amrex::Real vy = 6.0*y + 1.3*z - 2.9;
                amrex::Real vz = -4.6*z + 1.3*y + 3.8*x;
                
                vel_arr(i,j,k,2) = 4.0*x*x*z + 3.0*y*y - 2.3*z*z + 3.8*x*z + 9.4*x*y + 3.4*x - 2.9*y + 5.6;
                amrex::Real wx = 8.0*x*z + 3.8*z + 9.4*y + 3.4;
                amrex::Real wy = 6.0*y + 9.4*x -2.9;
                amrex::Real wz = 4.0*x*x -4.6*z + 3.8*x;

                sr_arr(i,j,k) = std::sqrt(2.0 * ux*ux + 2.0 * vy*vy + 2.0 * wz*wz + (uy+vx)*(uy+vx) + (vz+wy)*(vz+wy) + (wx+uz)*(wx+uz));
                
            }
        }
    }
    
    for(i = 1; i < n; ++i){
        for(j = 1; j < n; ++j){
            for(k = 1; k < n; ++k){
                sr = amr_wind::strainrate<amr_wind::StencilInterior>(i,j,k,idx,idy,idz,velocity.array());
                EXPECT_NEAR(sr_arr(i,j,k), sr, tol);
            }
        }
    }
    
}

TEST_F(StrainrateTest, ilo)
{

    constexpr double tol = 1.0e-12;
    
    const int m = 2;
    const int n = 4;
    amrex::Box bx{{0, 0, 0}, {m, n, n}};
    amrex::FArrayBox velocity(bx, AMREX_SPACEDIM);
    amrex::FArrayBox strainrate(bx, 1);

    velocity.setVal<amrex::RunOn::Host>(0.0);

    int i = 1, j = 2, k = 2;
    amrex::Real dx = 0.1, dy = 0.2, dz = 0.3 + 0.01*amrex::Random();
    amrex::Real idx = 1.0/dx, idy = 1.0/dy, idz = 1.0/dz;
    
    amrex::Real sr = amr_wind::strainrate<amr_wind::StencilILO>(i,j,k,idx,idy,idz,velocity.array());
    EXPECT_NEAR(sr, 0.0, tol);
    
    auto vel_arr = velocity.array();
    auto sr_arr = strainrate.array();
    
    // side stencil should be exact for linear polynomials
    for(i = 0; i <= m; ++i){
        for(j = 0; j <= n; ++j){
            for(k = 0; k <= n; ++k){
                
                amrex::Real x = 2.0;
                
                if(i==1) x += 0.5*dx;
                if(i==2) x += 1.5*dx;
                
                amrex::Real y = 2.5 + j*dy;
                amrex::Real z = 5.2 + k*dz;
                
                vel_arr(i,j,k,0) = 3.0*x - 2.6*y*z + 2.9*z*z + 7.8;
                amrex::Real ux =  3.0;
                amrex::Real uy = -2.6*z;
                amrex::Real uz = -2.6*y + 5.8*z;
                
                vel_arr(i,j,k,1) = 3.9*x + 0.6*y*y - 9.9*z + 3.8;
                amrex::Real vx = 3.9;
                amrex::Real vy = 1.2*y;
                amrex::Real vz = -9.9;
                
                vel_arr(i,j,k,2) = -5.6*x + 2.8*y + 4.5*z*y + 1.0;
                amrex::Real wx = -5.6;
                amrex::Real wy = 2.8 + 4.5*z;
                amrex::Real wz = 4.5*y;

                sr_arr(i,j,k) = std::sqrt(2.0 * ux*ux + 2.0 * vy*vy + 2.0 * wz*wz + (uy+vx)*(uy+vx) + (vz+wy)*(vz+wy) + (wx+uz)*(wx+uz));
                
            }
        }
    }
    
    i = 1;
    for(j = 1; j < n; ++j){
        for(k = 1; k < n; ++k){
            sr = amr_wind::strainrate<amr_wind::StencilILO>(i,j,k,idx,idy,idz,velocity.array());
            EXPECT_NEAR(sr_arr(i,j,k), sr, tol);
        }
    }
    
}

TEST_F(StrainrateTest, ihi)
{

    constexpr double tol = 1.0e-12;
    
    const int m = 2;
    const int n = 4;
    amrex::Box bx{{0, 0, 0}, {m, n, n}};
    amrex::FArrayBox velocity(bx, AMREX_SPACEDIM);
    amrex::FArrayBox strainrate(bx, 1);

    velocity.setVal<amrex::RunOn::Host>(0.0);

    int i = 1, j = 2, k = 2;
    amrex::Real dx = 0.1, dy = 0.2, dz = 0.3 + 0.01*amrex::Random();
    amrex::Real idx = 1.0/dx, idy = 1.0/dy, idz = 1.0/dz;
    
    amrex::Real sr = amr_wind::strainrate<amr_wind::StencilILO>(i,j,k,idx,idy,idz,velocity.array());
    EXPECT_NEAR(sr, 0.0, tol);
    
    auto vel_arr = velocity.array();
    auto sr_arr = strainrate.array();
    
    // side stencil should be exact for linear polynomials
    for(i = 0; i <= m; ++i){
        for(j = 0; j <= n; ++j){
            for(k = 0; k <= n; ++k){
                
                amrex::Real x = 2.0;
                
                if(i==1) x += dx;
                if(i==2) x += 1.5*dx;
                
                amrex::Real y = 2.5 + j*dy;
                amrex::Real z = 5.2 + k*dz;
                
                vel_arr(i,j,k,0) = 3.0*x - 2.6*y*z + 2.9*z*z + 7.8;
                amrex::Real ux =  3.0;
                amrex::Real uy = -2.6*z;
                amrex::Real uz = -2.6*y + 5.8*z;
                
                vel_arr(i,j,k,1) = 3.9*x + 0.6*y*y - 9.9*z + 3.8;
                amrex::Real vx = 3.9;
                amrex::Real vy = 1.2*y;
                amrex::Real vz = -9.9;
                
                vel_arr(i,j,k,2) = -5.6*x + 2.8*y + 4.5*z*y + 1.0;
                amrex::Real wx = -5.6;
                amrex::Real wy = 2.8 + 4.5*z;
                amrex::Real wz = 4.5*y;

                sr_arr(i,j,k) = std::sqrt(2.0 * ux*ux + 2.0 * vy*vy + 2.0 * wz*wz + (uy+vx)*(uy+vx) + (vz+wy)*(vz+wy) + (wx+uz)*(wx+uz));
                
            }
        }
    }
    
    i = 1;
    for(j = 1; j < n; ++j){
        for(k = 1; k < n; ++k){
            sr = amr_wind::strainrate<amr_wind::StencilIHI>(i,j,k,idx,idy,idz,velocity.array());
            EXPECT_NEAR(sr_arr(i,j,k), sr, tol);
        }
    }
    
}

TEST_F(StrainrateTest, jlo)
{

    constexpr double tol = 1.0e-12;
    
    const int m = 2;
    const int n = 4;
    amrex::Box bx{{0, 0, 0}, {n, m, n}};
    amrex::FArrayBox velocity(bx, AMREX_SPACEDIM);
    amrex::FArrayBox strainrate(bx, 1);

    velocity.setVal<amrex::RunOn::Host>(0.0);

    int i = 2, j = 1, k = 2;
    amrex::Real dx = 0.1, dy = 0.2, dz = 0.3 + 0.01*amrex::Random();
    amrex::Real idx = 1.0/dx, idy = 1.0/dy, idz = 1.0/dz;
    
    amrex::Real sr = amr_wind::strainrate<amr_wind::StencilJLO>(i,j,k,idx,idy,idz,velocity.array());
    EXPECT_NEAR(sr, 0.0, tol);
    
    auto vel_arr = velocity.array();
    auto sr_arr = strainrate.array();
    
    // side stencil should be exact for linear polynomials
    for(i = 0; i <= n; ++i){
        for(j = 0; j <= m; ++j){
            for(k = 0; k <= n; ++k){
                
                amrex::Real x = 2.0 + i*dx;
                amrex::Real y = 2.5;
                if(j==1) y += 0.5*dy;
                if(j==2) y += 1.5*dy;
                
                amrex::Real z = 5.2 + k*dz;
                
                vel_arr(i,j,k,0) = 3.0*x*x - 2.6*y + 2.9*z + 2.3;
                amrex::Real ux =  6.0*x;
                amrex::Real uy = -2.6;
                amrex::Real uz =  2.9;
                
                vel_arr(i,j,k,1) = 3.9*x*z + 0.6*y - 9.9*z + 3.8;
                amrex::Real vx = 3.9*z;
                amrex::Real vy = 0.6;
                amrex::Real vz = 3.9*x - 9.9;
                
                vel_arr(i,j,k,2) = -5.6*x + 2.8*y + 4.5*z + 2.3*z*z + 1.0;
                amrex::Real wx = -5.6;
                amrex::Real wy = 2.8;
                amrex::Real wz = 4.5 + 4.6*z;

                sr_arr(i,j,k) = std::sqrt(2.0 * ux*ux + 2.0 * vy*vy + 2.0 * wz*wz + (uy+vx)*(uy+vx) + (vz+wy)*(vz+wy) + (wx+uz)*(wx+uz));
                
            }
        }
    }
    
    j = 1;
    for(i = 1; i < n; ++i){
        for(k = 1; k < n; ++k){
            sr = amr_wind::strainrate<amr_wind::StencilJLO>(i,j,k,idx,idy,idz,velocity.array());
            EXPECT_NEAR(sr_arr(i,j,k), sr, tol);
        }
    }
    
}

TEST_F(StrainrateTest, jhi)
{

    constexpr double tol = 1.0e-12;
    
    const int m = 2;
    const int n = 4;
    amrex::Box bx{{0, 0, 0}, {n, m, n}};
    amrex::FArrayBox velocity(bx, AMREX_SPACEDIM);
    amrex::FArrayBox strainrate(bx, 1);

    velocity.setVal<amrex::RunOn::Host>(0.0);

    int i = 2, j = 1, k = 2;
    amrex::Real dx = 0.1, dy = 0.2, dz = 0.3 + 0.01*amrex::Random();
    amrex::Real idx = 1.0/dx, idy = 1.0/dy, idz = 1.0/dz;
    
    amrex::Real sr = amr_wind::strainrate<amr_wind::StencilJHI>(i,j,k,idx,idy,idz,velocity.array());
    EXPECT_NEAR(sr, 0.0, tol);
    
    auto vel_arr = velocity.array();
    auto sr_arr = strainrate.array();
    
    // side stencil should be exact for linear polynomials
    for(i = 0; i <= n; ++i){
        for(j = 0; j <= m; ++j){
            for(k = 0; k <= n; ++k){
                
                amrex::Real x = 2.0 + i*dx;
                amrex::Real y = 2.5;
                if(j==1) y += dy;
                if(j==2) y += 1.5*dy;
                amrex::Real z = 5.2 + k*dz;
                
                vel_arr(i,j,k,0) = 3.0*x*x - 2.6*y + 2.9*z + 2.3;
                amrex::Real ux =  6.0*x;
                amrex::Real uy = -2.6;
                amrex::Real uz =  2.9;
                
                vel_arr(i,j,k,1) = 3.9*x*z + 0.6*y - 9.9*z + 3.8;
                amrex::Real vx = 3.9*z;
                amrex::Real vy = 0.6;
                amrex::Real vz = 3.9*x - 9.9;
                
                vel_arr(i,j,k,2) = -5.6*x + 2.8*y + 4.5*z + 2.3*z*z + 1.0;
                amrex::Real wx = -5.6;
                amrex::Real wy = 2.8;
                amrex::Real wz = 4.5 + 4.6*z;

                sr_arr(i,j,k) = std::sqrt(2.0 * ux*ux + 2.0 * vy*vy + 2.0 * wz*wz + (uy+vx)*(uy+vx) + (vz+wy)*(vz+wy) + (wx+uz)*(wx+uz));
                
            }
        }
    }
    
    j = 1;
    for(i = 1; i < n; ++i){
        for(k = 1; k < n; ++k){
            sr = amr_wind::strainrate<amr_wind::StencilJHI>(i,j,k,idx,idy,idz,velocity.array());
            EXPECT_NEAR(sr_arr(i,j,k), sr, tol);
        }
    }
    
}

TEST_F(StrainrateTest, klo)
{

    constexpr double tol = 1.0e-12;
    
    const int m = 2;
    const int n = 4;
    amrex::Box bx{{0, 0, 0}, {n, n, m}};
    amrex::FArrayBox velocity(bx, AMREX_SPACEDIM);
    amrex::FArrayBox strainrate(bx, 1);

    velocity.setVal<amrex::RunOn::Host>(0.0);

    int i = 2, j = 2, k = 1;
    amrex::Real dx = 0.1, dy = 0.2, dz = 0.3 + 0.01*amrex::Random();
    amrex::Real idx = 1.0/dx, idy = 1.0/dy, idz = 1.0/dz;
    
    amrex::Real sr = amr_wind::strainrate<amr_wind::StencilKLO>(i,j,k,idx,idy,idz,velocity.array());
    EXPECT_NEAR(sr, 0.0, tol);
    
    auto vel_arr = velocity.array();
    auto sr_arr = strainrate.array();
    
    // side stencil should be exact for linear polynomials
    for(i = 0; i <= n; ++i){
        for(j = 0; j <= n; ++j){
            for(k = 0; k <= m; ++k){
                
                amrex::Real x = 2.0 + i*dx;
                amrex::Real y = 2.5 + j*dy;

                amrex::Real z = 5.2;
                if(k==1) z += 0.5*dz;
                if(k==2) z += 1.5*dz;
                
                vel_arr(i,j,k,0) = 3.0*x*x - 3.4*x*y + 2.6*y + 2.9*z + 2.3;
                amrex::Real ux =  6.0*x - 3.4*y;
                amrex::Real uy =  -3.4*x + 2.6;
                amrex::Real uz =  2.9;
                
                vel_arr(i,j,k,1) = 3.9*x + 0.6*y*y - 9.9*z + 3.8;
                amrex::Real vx = 3.9;
                amrex::Real vy = 1.2*y;
                amrex::Real vz = -9.9;
                
                vel_arr(i,j,k,2) = -5.6*x*y + 2.8*y + 4.5*z + 1.0;
                amrex::Real wx = -5.6*y;
                amrex::Real wy = -5.6*x + 2.8;
                amrex::Real wz = 4.5;

                sr_arr(i,j,k) = std::sqrt(2.0 * ux*ux + 2.0 * vy*vy + 2.0 * wz*wz + (uy+vx)*(uy+vx) + (vz+wy)*(vz+wy) + (wx+uz)*(wx+uz));
                
            }
        }
    }
    
    k = 1;
    for(i = 1; i < n; ++i){
        for(j = 1; j < n; ++j){
            sr = amr_wind::strainrate<amr_wind::StencilKLO>(i,j,k,idx,idy,idz,velocity.array());
            EXPECT_NEAR(sr_arr(i,j,k), sr, tol);
        }
    }
    
}

TEST_F(StrainrateTest, khi)
{

    constexpr double tol = 1.0e-12;
    
    const int m = 2;
    const int n = 4;
    amrex::Box bx{{0, 0, 0}, {n, n, m}};
    amrex::FArrayBox velocity(bx, AMREX_SPACEDIM);
    amrex::FArrayBox strainrate(bx, 1);

    velocity.setVal<amrex::RunOn::Host>(0.0);

    int i = 2, j = 2, k = 1;
    amrex::Real dx = 0.1, dy = 0.2, dz = 0.3 + 0.01*amrex::Random();
    amrex::Real idx = 1.0/dx, idy = 1.0/dy, idz = 1.0/dz;
    
    amrex::Real sr = amr_wind::strainrate<amr_wind::StencilKHI>(i,j,k,idx,idy,idz,velocity.array());
    EXPECT_NEAR(sr, 0.0, tol);
    
    auto vel_arr = velocity.array();
    auto sr_arr = strainrate.array();
    
    // side stencil should be exact for linear polynomials
    for(i = 0; i <= n; ++i){
        for(j = 0; j <= n; ++j){
            for(k = 0; k <= m; ++k){
                
                amrex::Real x = 2.0 + i*dx;
                amrex::Real y = 2.5 + j*dy;
                amrex::Real z = 5.2;
                if(k==1) z += dz;
                if(k==2) z += 1.5*dz;
                
                vel_arr(i,j,k,0) = 3.0*x*x - 3.4*x*y + 2.6*y + 2.9*z + 2.3;
                amrex::Real ux =  6.0*x - 3.4*y;
                amrex::Real uy =  -3.4*x + 2.6;
                amrex::Real uz =  2.9;
                
                vel_arr(i,j,k,1) = 3.9*x + 0.6*y*y - 9.9*z + 3.8;
                amrex::Real vx = 3.9;
                amrex::Real vy = 1.2*y;
                amrex::Real vz = -9.9;
                
                vel_arr(i,j,k,2) = -5.6*x*y + 2.8*y + 4.5*z + 1.0;
                amrex::Real wx = -5.6*y;
                amrex::Real wy = -5.6*x + 2.8;
                amrex::Real wz = 4.5;

                sr_arr(i,j,k) = std::sqrt(2.0 * ux*ux + 2.0 * vy*vy + 2.0 * wz*wz + (uy+vx)*(uy+vx) + (vz+wy)*(vz+wy) + (wx+uz)*(wx+uz));
                
            }
        }
    }
    
    k = 1;
    for(i = 1; i < n; ++i){
        for(j = 1; j < n; ++j){
            sr = amr_wind::strainrate<amr_wind::StencilKHI>(i,j,k,idx,idy,idz,velocity.array());
            EXPECT_NEAR(sr_arr(i,j,k), sr, tol);
        }
    }
    
}

} // namespace amr_wind_tests
