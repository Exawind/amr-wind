#ifndef DERIVE_K_H_
#define DERIVE_K_H_

#include "Field.H"
#include <AMReX_FArrayBox.H>
#include <cmath>

namespace amr_wind {

struct StencilInterior
{
    // First derivatives
    static constexpr amrex::Real c00 =  0.5;
    static constexpr amrex::Real c01 =  0.0;
    static constexpr amrex::Real c02 = -0.5;
    static constexpr amrex::Real c10 =  c00;
    static constexpr amrex::Real c11 =  c01;
    static constexpr amrex::Real c12 =  c02;
    static constexpr amrex::Real c20 =  c00;
    static constexpr amrex::Real c21 =  c01;
    static constexpr amrex::Real c22 =  c02;

    // Second derivatives
    static constexpr amrex::Real s00 =  1.0;
    static constexpr amrex::Real s01 = -2.0;
    static constexpr amrex::Real s02 =  1.0;
    static constexpr amrex::Real s10 =  s00;
    static constexpr amrex::Real s11 =  s01;
    static constexpr amrex::Real s12 =  s02;
    static constexpr amrex::Real s20 =  s00;
    static constexpr amrex::Real s21 =  s01;
    static constexpr amrex::Real s22 =  s02;

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



template <typename Stencil>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void compute_gradient_scalar(
    int i,
    int j,
    int k,
    amrex::Real idx,
    amrex::Real idy,
    amrex::Real idz,
    amrex::Array4<amrex::Real const> const& phi,
    amrex::Array4<amrex::Real> const& gradphi) noexcept
{
    using namespace amrex;

    Real cp1, c, cm1;

    cp1 = Stencil::c00;
    c   = Stencil::c01;
    cm1 = Stencil::c02;

    gradphi(i,j,k,0) = (cp1*phi(i+1,j,k) + c*phi(i,j,k) + cm1*phi(i-1,j,k)) * idx;
    cp1 = Stencil::c10;
    c   = Stencil::c11;
    cm1 = Stencil::c12;
    gradphi(i,j,k,1) = (cp1*phi(i,j+1,k) + c*phi(i,j,k) + cm1*phi(i,j-1,k)) * idy;
    cp1 = Stencil::c20;
    c   = Stencil::c21;
    cm1 = Stencil::c22;
    gradphi(i,j,k,2) = (cp1*phi(i,j,k+1) + c*phi(i,j,k) + cm1*phi(i,j,k-1)) * idz;

}

template <typename Stencil>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real curvature(
    int i,
    int j,
    int k,
    amrex::Real idx,
    amrex::Real idy,
    amrex::Real idz,
    amrex::Array4<amrex::Real const> const& phi,
    amrex::Array4<amrex::Real const> const& gradphi) noexcept
{
    using namespace amrex;

    Real cp1, c, cm1;

    cp1 = Stencil::c00;
    c   = Stencil::c01;
    cm1 = Stencil::c02;
    const Real phix= gradphi(i,j,k,0);
    const Real phixx = (cp1*gradphi(i+1,j,k,0) + c*gradphi(i,j,k,0) + cm1*gradphi(i-1,j,k,0)) * idx;
    cp1 = Stencil::c10;
    c   = Stencil::c11;
    cm1 = Stencil::c12;
    const Real phiy= gradphi(i,j,k,1);
    const Real phixy = (cp1*gradphi(i,j+1,k,0) + c*gradphi(i,j,k,0) + cm1*gradphi(i,j-1,k,0)) * idy;
    const Real phiyy = (cp1*gradphi(i,j+1,k,1) + c*gradphi(i,j,k,1) + cm1*gradphi(i,j-1,k,1)) * idy;
    cp1 = Stencil::c20;
    c   = Stencil::c21;
    cm1 = Stencil::c22;
    const Real phiz= gradphi(i,j,k,2);
    const Real phixz = (cp1*gradphi(i,j,k+1,0) + c*gradphi(i,j,k,0) + cm1*gradphi(i,j,k-1,0)) * idz;
    const Real phiyz = (cp1*gradphi(i,j,k+1,1) + c*gradphi(i,j,k,1) + cm1*gradphi(i,j,k-1,1)) * idz;
    const Real phizz = (cp1*gradphi(i,j,k+1,2) + c*gradphi(i,j,k,2) + cm1*gradphi(i,j,k-1,2)) * idz;
   
    const Real abs_gradphi=std::sqrt(gradphi(i,j,k,0)*gradphi(i,j,k,0)+
                          gradphi(i,j,k,1)*gradphi(i,j,k,1)+
                          gradphi(i,j,k,2)*gradphi(i,j,k,2));
    return (phix*phix*phiyy - 2.*phix*phiy*phixy + phiy*phiy*phixx + phix*phix*phizz
            -2.*phix*phiz*phixz + phiz*phiz*phixx + phiy*phiy*phizz - 2*phiy*phiz*phiyz 
            + phiz*phiz*phiyy)/(abs_gradphi*abs_gradphi*abs_gradphi);
}



template <typename Stencil>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real strainrate(
    int i,
    int j,
    int k,
    amrex::Real idx,
    amrex::Real idy,
    amrex::Real idz,
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

template<typename FType>
void compute_strainrate(FType& field, const Field& velocity);


} // namespace amr_wind
#endif
