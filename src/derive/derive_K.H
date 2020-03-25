#ifndef DERIVE_K_H_
#define DERIVE_K_H_

#include <AMReX_FArrayBox.H>
#include <cmath>

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real incflo_strainrate (int i, int j, int k,
                               amrex::Real idx, amrex::Real idy, amrex::Real idz,
                               amrex::Array4<amrex::Real const> const& vel) noexcept
{
    using namespace amrex;
    
    const Real ux = 0.5 * (vel(i+1,j,k,0) - vel(i-1,j,k,0)) * idx;
    const Real vx = 0.5 * (vel(i+1,j,k,1) - vel(i-1,j,k,1)) * idx;
    const Real wx = 0.5 * (vel(i+1,j,k,2) - vel(i-1,j,k,2)) * idx;

    const Real uy = 0.5 * (vel(i,j+1,k,0) - vel(i,j-1,k,0)) * idy;
    const Real vy = 0.5 * (vel(i,j+1,k,1) - vel(i,j-1,k,1)) * idy;
    const Real wy = 0.5 * (vel(i,j+1,k,2) - vel(i,j-1,k,2)) * idy;

    const Real uz = 0.5 * (vel(i,j,k+1,0) - vel(i,j,k-1,0)) * idz;
    const Real vz = 0.5 * (vel(i,j,k+1,1) - vel(i,j,k-1,1)) * idz;
    const Real wz = 0.5 * (vel(i,j,k+1,2) - vel(i,j,k-1,2)) * idz;

    return std::sqrt(2.0 * ux*ux + 2.0 * vy*vy + 2.0 * wz*wz
                     + (uy+vx)*(uy+vx) + (vz+wy)*(vz+wy) + (wx+uz)*(wx+uz));
}

// i = 0 first cell
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real incflo_strainrate_ilo (int i, int j, int k,
                               amrex::Real idx, amrex::Real idy, amrex::Real idz,
                               amrex::Array4<amrex::Real const> const& vel) noexcept
{
    using namespace amrex;

    // average velocity to face i+1/2 and take derivative about cell center between interior face and boundary face i-1
    const Real ux = (0.5 * (vel(i+1,j,k,0) + vel(i,j,k,0)) - vel(i-1,j,k,0)) * idx;
    const Real vx = (0.5 * (vel(i+1,j,k,1) + vel(i,j,k,1)) - vel(i-1,j,k,1)) * idx;
    const Real wx = (0.5 * (vel(i+1,j,k,2) + vel(i,j,k,2)) - vel(i-1,j,k,2)) * idx;
    
    const Real uy = 0.5 * (vel(i,j+1,k,0) - vel(i,j-1,k,0)) * idy;
    const Real vy = 0.5 * (vel(i,j+1,k,1) - vel(i,j-1,k,1)) * idy;
    const Real wy = 0.5 * (vel(i,j+1,k,2) - vel(i,j-1,k,2)) * idy;

    const Real uz = 0.5 * (vel(i,j,k+1,0) - vel(i,j,k-1,0)) * idz;
    const Real vz = 0.5 * (vel(i,j,k+1,1) - vel(i,j,k-1,1)) * idz;
    const Real wz = 0.5 * (vel(i,j,k+1,2) - vel(i,j,k-1,2)) * idz;

    return std::sqrt(2.0 * ux*ux + 2.0 * vy*vy + 2.0 * wz*wz
                     + (uy+vx)*(uy+vx) + (vz+wy)*(vz+wy) + (wx+uz)*(wx+uz));
}

// j = 0 first cell
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real incflo_strainrate_jlo (int i, int j, int k,
                               amrex::Real idx, amrex::Real idy, amrex::Real idz,
                               amrex::Array4<amrex::Real const> const& vel) noexcept
{
    using namespace amrex;

    const Real ux = 0.5 * (vel(i+1,j,k,0) - vel(i-1,j,k,0)) * idx;
    const Real vx = 0.5 * (vel(i+1,j,k,1) - vel(i-1,j,k,1)) * idx;
    const Real wx = 0.5 * (vel(i+1,j,k,2) - vel(i-1,j,k,2)) * idx;

    // average velocity to face k+1/2 and take derivative about cell center between interior face and boundary face k-1
    const Real uy = (0.5 * (vel(i,j+1,k,0) + vel(i,j,k,0)) - vel(i,j-1,k,0)) * idy;
    const Real vy = (0.5 * (vel(i,j+1,k,1) + vel(i,j,k,1)) - vel(i,j-1,k,1)) * idy;
    const Real wy = (0.5 * (vel(i,j+1,k,2) + vel(i,j,k,2)) - vel(i,j-1,k,2)) * idy;
   
    const Real uz = 0.5 * (vel(i,j,k+1,0) - vel(i,j,k-1,0)) * idz;
    const Real vz = 0.5 * (vel(i,j,k+1,1) - vel(i,j,k-1,1)) * idz;
    const Real wz = 0.5 * (vel(i,j,k+1,2) - vel(i,j,k-1,2)) * idz;
    
    return std::sqrt(2.0 * ux*ux + 2.0 * vy*vy + 2.0 * wz*wz
                     + (uy+vx)*(uy+vx) + (vz+wy)*(vz+wy) + (wx+uz)*(wx+uz));
}

// k = 0 first cell
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real incflo_strainrate_klo (int i, int j, int k,
                               amrex::Real idx, amrex::Real idy, amrex::Real idz,
                               amrex::Array4<amrex::Real const> const& vel) noexcept
{
    using namespace amrex;

    const Real ux = 0.5 * (vel(i+1,j,k,0) - vel(i-1,j,k,0)) * idx;
    const Real vx = 0.5 * (vel(i+1,j,k,1) - vel(i-1,j,k,1)) * idx;
    const Real wx = 0.5 * (vel(i+1,j,k,2) - vel(i-1,j,k,2)) * idx;

    const Real uy = 0.5 * (vel(i,j+1,k,0) - vel(i,j-1,k,0)) * idy;
    const Real vy = 0.5 * (vel(i,j+1,k,1) - vel(i,j-1,k,1)) * idy;
    const Real wy = 0.5 * (vel(i,j+1,k,2) - vel(i,j-1,k,2)) * idy;

    // average velocity to face k+1/2 and take derivative about cell center between interior face and boundary face k-1
    const Real uz = (0.5 * (vel(i,j,k+1,0) + vel(i,j,k,0)) - vel(i,j,k-1,0)) * idz;
    const Real vz = (0.5 * (vel(i,j,k+1,1) + vel(i,j,k,1)) - vel(i,j,k-1,1)) * idz;
    const Real wz = (0.5 * (vel(i,j,k+1,2) + vel(i,j,k,2)) - vel(i,j,k-1,2)) * idz;

    return std::sqrt(2.0 * ux*ux + 2.0 * vy*vy + 2.0 * wz*wz
                     + (uy+vx)*(uy+vx) + (vz+wy)*(vz+wy) + (wx+uz)*(wx+uz));
}

// i = ncell-1
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real incflo_strainrate_ihi (int i, int j, int k,
                               amrex::Real idx, amrex::Real idy, amrex::Real idz,
                               amrex::Array4<amrex::Real const> const& vel) noexcept
{
    using namespace amrex;

    // average velocity to face i-1/2 and take derivative about cell center between interior face and boundary face i+1
    const Real ux = (vel(i+1,j,k,0) - 0.5 * (vel(i,j,k,0) + vel(i-1,j,k,0))) * idx;
    const Real vx = (vel(i+1,j,k,1) - 0.5 * (vel(i,j,k,1) + vel(i-1,j,k,1))) * idx;
    const Real wx = (vel(i+1,j,k,2) - 0.5 * (vel(i,j,k,2) + vel(i-1,j,k,2))) * idx;

    const Real uy = 0.5 * (vel(i,j+1,k,0) - vel(i,j-1,k,0)) * idy;
    const Real vy = 0.5 * (vel(i,j+1,k,1) - vel(i,j-1,k,1)) * idy;
    const Real wy = 0.5 * (vel(i,j+1,k,2) - vel(i,j-1,k,2)) * idy;

    const Real uz = 0.5 * (vel(i,j,k+1,0) - vel(i,j,k-1,0)) * idz;
    const Real vz = 0.5 * (vel(i,j,k+1,1) - vel(i,j,k-1,1)) * idz;
    const Real wz = 0.5 * (vel(i,j,k+1,2) - vel(i,j,k-1,2)) * idz;

    return std::sqrt(2.0 * ux*ux + 2.0 * vy*vy + 2.0 * wz*wz
                     + (uy+vx)*(uy+vx) + (vz+wy)*(vz+wy) + (wx+uz)*(wx+uz));
}

// j = ncell-1
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real incflo_strainrate_jhi (int i, int j, int k,
                               amrex::Real idx, amrex::Real idy, amrex::Real idz,
                               amrex::Array4<amrex::Real const> const& vel) noexcept
{
    using namespace amrex;

    const Real ux = 0.5 * (vel(i+1,j,k,0) - vel(i-1,j,k,0)) * idx;
    const Real vx = 0.5 * (vel(i+1,j,k,1) - vel(i-1,j,k,1)) * idx;
    const Real wx = 0.5 * (vel(i+1,j,k,2) - vel(i-1,j,k,2)) * idx;

    // average velocity to face j-1/2 and take derivative about cell center between interior face and boundary face j+1
    const Real uy = (vel(i,j+1,k,0) - 0.5 * (vel(i,j,k,0) + vel(i,j-1,k,0))) * idy;
    const Real vy = (vel(i,j+1,k,1) - 0.5 * (vel(i,j,k,1) + vel(i,j-1,k,1))) * idy;
    const Real wy = (vel(i,j+1,k,2) - 0.5 * (vel(i,j,k,2) + vel(i,j-1,k,2))) * idy;

    const Real uz = 0.5 * (vel(i,j,k+1,0) - vel(i,j,k-1,0)) * idz;
    const Real vz = 0.5 * (vel(i,j,k+1,1) - vel(i,j,k-1,1)) * idz;
    const Real wz = 0.5 * (vel(i,j,k+1,2) - vel(i,j,k-1,2)) * idz;

    return std::sqrt(2.0 * ux*ux + 2.0 * vy*vy + 2.0 * wz*wz
                     + (uy+vx)*(uy+vx) + (vz+wy)*(vz+wy) + (wx+uz)*(wx+uz));
}

// k = ncell-1
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real incflo_strainrate_khi (int i, int j, int k,
                               amrex::Real idx, amrex::Real idy, amrex::Real idz,
                               amrex::Array4<amrex::Real const> const& vel) noexcept
{
    using namespace amrex;

    const Real ux = 0.5 * (vel(i+1,j,k,0) - vel(i-1,j,k,0)) * idx;
    const Real vx = 0.5 * (vel(i+1,j,k,1) - vel(i-1,j,k,1)) * idx;
    const Real wx = 0.5 * (vel(i+1,j,k,2) - vel(i-1,j,k,2)) * idx;

    const Real uy = 0.5 * (vel(i,j+1,k,0) - vel(i,j-1,k,0)) * idy;
    const Real vy = 0.5 * (vel(i,j+1,k,1) - vel(i,j-1,k,1)) * idy;
    const Real wy = 0.5 * (vel(i,j+1,k,2) - vel(i,j-1,k,2)) * idy;

    // average velocity to face k-1/2 and take derivative about cell center between interior face and boundary face k+1
    const Real uz = (vel(i,j,k+1,0) - 0.5 * (vel(i,j,k,0) + vel(i,j,k-1,0))) * idz;
    const Real vz = (vel(i,j,k+1,1) - 0.5 * (vel(i,j,k,1) + vel(i,j,k-1,1))) * idz;
    const Real wz = (vel(i,j,k+1,2) - 0.5 * (vel(i,j,k,2) + vel(i,j,k-1,2))) * idz;

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
