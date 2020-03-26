#ifndef DERIVE_K_H_
#define DERIVE_K_H_

#include <AMReX_FArrayBox.H>
#include <cmath>

struct StencilInterior
{
    static constexpr amrex::Real c00 =  0.5;
    static constexpr amrex::Real c01 =  0.0;
    static constexpr amrex::Real c02 = -0.5;
    static constexpr amrex::Real c10 =  c00;
    static constexpr amrex::Real c11 =  c01;
    static constexpr amrex::Real c12 =  c02;
    static constexpr amrex::Real c20 =  c00;
    static constexpr amrex::Real c21 =  c01;
    static constexpr amrex::Real c22 =  c02;
};

struct StencilILO
{
    static constexpr amrex::Real c00 =  0.5;
    static constexpr amrex::Real c01 =  0.5;
    static constexpr amrex::Real c02 = -1.0;
    static constexpr amrex::Real c10 =  0.5;
    static constexpr amrex::Real c11 =  0.0;
    static constexpr amrex::Real c12 = -0.5;
    static constexpr amrex::Real c20 =  c10;
    static constexpr amrex::Real c21 =  c11;
    static constexpr amrex::Real c22 =  c12;
};

struct StencilJLO
{
    static constexpr amrex::Real c00 =  0.5;
    static constexpr amrex::Real c01 =  0.0;
    static constexpr amrex::Real c02 = -0.5;
    static constexpr amrex::Real c10 =  0.5;
    static constexpr amrex::Real c11 =  0.5;
    static constexpr amrex::Real c12 = -1.0;
    static constexpr amrex::Real c20 =  c00;
    static constexpr amrex::Real c21 =  c01;
    static constexpr amrex::Real c22 =  c02;
};

struct StencilKLO
{
    static constexpr amrex::Real c00 =  0.5;
    static constexpr amrex::Real c01 =  0.0;
    static constexpr amrex::Real c02 = -0.5;
    static constexpr amrex::Real c10 =  c00;
    static constexpr amrex::Real c11 =  c01;
    static constexpr amrex::Real c12 =  c02;
    static constexpr amrex::Real c20 =  0.5;
    static constexpr amrex::Real c21 =  0.5;
    static constexpr amrex::Real c22 = -1.0;
};

struct StencilIHI
{
    static constexpr amrex::Real c00 =  1.0;
    static constexpr amrex::Real c01 = -0.5;
    static constexpr amrex::Real c02 = -0.5;
    static constexpr amrex::Real c10 =  0.5;
    static constexpr amrex::Real c11 =  0.0;
    static constexpr amrex::Real c12 = -0.5;
    static constexpr amrex::Real c20 =  c10;
    static constexpr amrex::Real c21 =  c11;
    static constexpr amrex::Real c22 =  c12;
};

struct StencilJHI
{
    static constexpr amrex::Real c00 =  0.5;
    static constexpr amrex::Real c01 =  0.0;
    static constexpr amrex::Real c02 = -0.5;
    static constexpr amrex::Real c10 =  1.0;
    static constexpr amrex::Real c11 = -0.5;
    static constexpr amrex::Real c12 = -0.5;
    static constexpr amrex::Real c20 =  c00;
    static constexpr amrex::Real c21 =  c01;
    static constexpr amrex::Real c22 =  c02;
};

struct StencilKHI
{
    static constexpr amrex::Real c00 =  0.5;
    static constexpr amrex::Real c01 =  0.0;
    static constexpr amrex::Real c02 = -0.5;
    static constexpr amrex::Real c10 =  c00;
    static constexpr amrex::Real c11 =  c01;
    static constexpr amrex::Real c12 =  c02;
    static constexpr amrex::Real c20 =  1.0;
    static constexpr amrex::Real c21 = -0.5;
    static constexpr amrex::Real c22 = -0.5;
};

template<typename Stencil>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real incflo_strainrate (int i, int j, int k,
                               amrex::Real idx, amrex::Real idy, amrex::Real idz,
                               amrex::Array4<amrex::Real const> const& vel) noexcept
{
    using namespace amrex;
    
    Real cp1, c, cm1;
    
    cp1 = Stencil::c00;
    c   = Stencil::c01;
    cm1 = Stencil::c02;

    const Real ux = (cp1*vel(i+1,j,k,0) + c*vel(i,j,k,0) + cm1*vel(i-1,j,k,0)) * idx;
    const Real vx = (cp1*vel(i+1,j,k,1) + c*vel(i,j,k,1) + cm1*vel(i-1,j,k,1)) * idx;
    const Real wx = (cp1*vel(i+1,j,k,2) + c*vel(i,j,k,2) + cm1*vel(i-1,j,k,2)) * idx;

    cp1 = Stencil::c10;
    c   = Stencil::c11;
    cm1 = Stencil::c12;
    
    const Real uy = (cp1*vel(i,j+1,k,0) + c*vel(i,j,k,0) + cm1*vel(i,j-1,k,0)) * idy;
    const Real vy = (cp1*vel(i,j+1,k,1) + c*vel(i,j,k,1) + cm1*vel(i,j-1,k,1)) * idy;
    const Real wy = (cp1*vel(i,j+1,k,2) + c*vel(i,j,k,2) + cm1*vel(i,j-1,k,2)) * idy;

    cp1 = Stencil::c20;
    c   = Stencil::c21;
    cm1 = Stencil::c22;
    
    const Real uz = (cp1*vel(i,j,k+1,0) + c*vel(i,j,k,0) + cm1*vel(i,j,k-1,0)) * idz;
    const Real vz = (cp1*vel(i,j,k+1,1) + c*vel(i,j,k,1) + cm1*vel(i,j,k-1,1)) * idz;
    const Real wz = (cp1*vel(i,j,k+1,2) + c*vel(i,j,k,2) + cm1*vel(i,j,k-1,2)) * idz;

    return std::sqrt(2.0 * ux*ux + 2.0 * vy*vy + 2.0 * wz*wz
                     + (uy+vx)*(uy+vx) + (vz+wy)*(vz+wy) + (wx+uz)*(wx+uz));
}

#ifdef AMREX_USE_EB
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real incflo_strainrate_eb (int i, int j, int k,
                                  amrex::Real idx, amrex::Real idy, amrex::Real idz,
                                  amrex::Array4<amrex::Real const> const& vel,
                                  amrex::EBCellFlag flag) noexcept
{
    using namespace amrex;

    if (flag.isRegular()) {
        return incflo_strainrate(i,j,k,idx,idy,idz,vel);
    } else if (flag.isSingleValued()) {
        constexpr Real c0 = -1.5;
        constexpr Real c1 =  2.0;
        constexpr Real c2 = -0.5;
        Real ux, vx, wx, uy, vy, wy, uz, vz, wz;
        if (!flag.isConnected( 1,0,0))
        {
            // Covered cell to the right, go fish left
            ux = - (c0 * vel(i  ,j,k,0)
                  + c1 * vel(i-1,j,k,0)
                  + c2 * vel(i-2,j,k,0)) * idx;
            vx = - (c0 * vel(i  ,j,k,1)
                  + c1 * vel(i-1,j,k,1)
                  + c2 * vel(i-2,j,k,1)) * idx;
            wx = - (c0 * vel(i  ,j,k,2)
                  + c1 * vel(i-1,j,k,2)
                  + c2 * vel(i-2,j,k,2)) * idx;
        }
        else if (!flag.isConnected(-1,0,0))
        {
            // Covered cell to the left, go fish right
            ux = (c0 * vel(i  ,j,k,0)
                + c1 * vel(i+1,j,k,0)
                + c2 * vel(i+2,j,k,0)) * idx;
            vx = (c0 * vel(i  ,j,k,1)
                + c1 * vel(i+1,j,k,1)
                + c2 * vel(i+2,j,k,1)) * idx;
            wx = (c0 * vel(i  ,j,k,2)
                + c1 * vel(i+1,j,k,2)
                + c2 * vel(i+2,j,k,2)) * idx;
        }
        else
        {
            // No covered cells right or left, use standard stencil
            ux = 0.5 * (vel(i+1,j,k,0) - vel(i-1,j,k,0)) * idx;
            vx = 0.5 * (vel(i+1,j,k,1) - vel(i-1,j,k,1)) * idx;
            wx = 0.5 * (vel(i+1,j,k,2) - vel(i-1,j,k,2)) * idx;
        }

        // Do the same in y-direction
        if (!flag.isConnected(0, 1,0))
        {
            uy = - (c0 * vel(i,j  ,k,0)
                  + c1 * vel(i,j-1,k,0)
                  + c2 * vel(i,j-2,k,0)) * idy;
            vy = - (c0 * vel(i,j  ,k,1)
                  + c1 * vel(i,j-1,k,1)
                  + c2 * vel(i,j-2,k,1)) * idy;
            wy = - (c0 * vel(i,j  ,k,2)
                  + c1 * vel(i,j-1,k,2)
                  + c2 * vel(i,j-2,k,2)) * idy;
        }
        else if (!flag.isConnected(0,-1,0))
        {
            uy = (c0 * vel(i,j  ,k,0)
                + c1 * vel(i,j+1,k,0)
                + c2 * vel(i,j+2,k,0)) * idy;
            vy = (c0 * vel(i,j  ,k,1)
                + c1 * vel(i,j+1,k,1)
                + c2 * vel(i,j+2,k,1)) * idy;
            wy = (c0 * vel(i,j  ,k,2)
                + c1 * vel(i,j+1,k,2)
                + c2 * vel(i,j+2,k,2)) * idy;
        }
        else
        {
            uy = 0.5 * (vel(i,j+1,k,0) - vel(i,j-1,k,0)) * idy;
            vy = 0.5 * (vel(i,j+1,k,1) - vel(i,j-1,k,1)) * idy;
            wy = 0.5 * (vel(i,j+1,k,2) - vel(i,j-1,k,2)) * idy;
        }

        // Do the same in z-direction
        if (!flag.isConnected(0,0, 1))
        {
            uz = - (c0 * vel(i,j,k  ,0)
                  + c1 * vel(i,j,k-1,0)
                  + c2 * vel(i,j,k-2,0)) * idz;
            vz = - (c0 * vel(i,j,k  ,1)
                  + c1 * vel(i,j,k-1,1)
                  + c2 * vel(i,j,k-2,1)) * idz;
            wz = - (c0 * vel(i,j,k  ,2)
                  + c1 * vel(i,j,k-1,2)
                  + c2 * vel(i,j,k-2,2)) * idz;
        }
        else if (!flag.isConnected(0,0,-1))
        {
            uz = (c0 * vel(i,j,k  ,0)
                + c1 * vel(i,j,k+1,0)
                + c2 * vel(i,j,k+2,0)) * idz;
            vz = (c0 * vel(i,j,k  ,1)
                + c1 * vel(i,j,k+1,1)
                + c2 * vel(i,j,k+2,1)) * idz;
            wz = (c0 * vel(i,j,k  ,2)
                + c1 * vel(i,j,k+1,2)
                + c2 * vel(i,j,k+2,2)) * idz;
        }
        else
        {
            uz = 0.5 * (vel(i,j,k+1,0) - vel(i,j,k-1,0)) * idz;
            vz = 0.5 * (vel(i,j,k+1,1) - vel(i,j,k-1,1)) * idz;
            wz = 0.5 * (vel(i,j,k+1,2) - vel(i,j,k-1,2)) * idz;
        }

        return std::sqrt(2.0 * ux*ux + 2.0 * vy*vy + 2.0 * wz*wz
                         + (uy+vx)*(uy+vx) + (vz+wy)*(vz+wy) + (wx+uz)*(wx+uz));
    } else {
        return 0.0;
    }
}
#endif

#endif
