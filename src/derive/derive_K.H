#ifndef DERIVE_K_H_
#define DERIVE_K_H_

#include "Field.H"
#include <AMReX_FArrayBox.H>
#include <cmath>

namespace amr_wind {

struct StencilInterior
{
    // First derivatives
    // fx= (ax*f_{i+1,j,k}+bx*f_{i,j,k}+cx*f_{i-1,j,k})/dx
    // fy= (ay*f_{i,j+1,k}+by*f_{i,j,k}+cy*f_{i,j-1,k})/dy
    // fz= (az*f_{i,j,k+1}+bz*f_{i,j,k}+cz*f_{i,j,k-1})/dz
    static constexpr amrex::Real c00 =  0.5; // ax 
    static constexpr amrex::Real c01 =  0.0; // bx
    static constexpr amrex::Real c02 = -0.5; // cx
    static constexpr amrex::Real c10 =  c00; // ay
    static constexpr amrex::Real c11 =  c01; // by
    static constexpr amrex::Real c12 =  c02; // cy
    static constexpr amrex::Real c20 =  c00; // az
    static constexpr amrex::Real c21 =  c01; // bz
    static constexpr amrex::Real c22 =  c02; // cz

    // Second derivatives
    // fxx= (axx*f_{i+1,j,k}+bxx*f_{i,j,k}+cxx*f_{i-1,j,k})/dx^2
    // fyy= (ayy*f_{i,j+1,k}+byy*f_{i,j,k}+cyy*f_{i,j-1,k})/dy^2
    // fzz= (azz*f_{i,j,k+1}+bzz*f_{i,j,k}+czz*f_{i,j,k-1})/dz^2
    static constexpr amrex::Real s00 =  1.0; // axx 
    static constexpr amrex::Real s01 = -2.0; // bxx
    static constexpr amrex::Real s02 =  1.0; // cxx
    static constexpr amrex::Real s10 =  s00; // ayy
    static constexpr amrex::Real s11 =  s01; // byy
    static constexpr amrex::Real s12 =  s02; // cyy
    static constexpr amrex::Real s20 =  s00; // azz
    static constexpr amrex::Real s21 =  s01; // bzz
    static constexpr amrex::Real s22 =  s02; // czz

    // Mixed derivatives
    // fxy = (axy*f_{i+1,j+1}+bxy*f_{i-1,j+1}+cxy*f_{i+1,j-1}+dxy*f_{i+1,j-1}+exy*f_{i-1,j-1}
    // fxz = (axz*f_{i+1,j+1}+bxz*f_{i-1,j+1}+cxz*f_{i+1,j-1}+dxz*f_{i+1,j-1}+exz*f_{i-1,j-1}
    // fyz = (ayz*f_{i+1,j+1}+byz*f_{i-1,j+1}+cyz*f_{i+1,j-1}+dyz*f_{i+1,j-1}+eyz*f_{i-1,j-1}
    static constexpr amrex::Real m00  =   0.25; // axy 
    static constexpr amrex::Real m01  =  -0.25; // bxy
    static constexpr amrex::Real m02  =   0.;   // cxy
    static constexpr amrex::Real m03  =  -0.25; // dxy
    static constexpr amrex::Real m04  =  0.25; // 3xy
    static constexpr amrex::Real m10  =  m00;   // axz 
    static constexpr amrex::Real m11  =  m01;   // bxz
    static constexpr amrex::Real m12  =  m02;   // cxz
    static constexpr amrex::Real m13  =  m03;   // dyz
    static constexpr amrex::Real m14  =  m04;   // eyz
    static constexpr amrex::Real m20  =  m00;   // ayz 
    static constexpr amrex::Real m21  =  m01;   // byz
    static constexpr amrex::Real m22  =  m02;   // cyz
    static constexpr amrex::Real m23  =  m03;   // dyz
    static constexpr amrex::Real m24  =  m04;   // eyz
};

struct StencilILO
{
    static constexpr amrex::Real c00 =  1./3.;
    static constexpr amrex::Real c01 =  1.;
    static constexpr amrex::Real c02 = -4./3.;
    static constexpr amrex::Real c10 =  0.5;
    static constexpr amrex::Real c11 =  0.0;
    static constexpr amrex::Real c12 = -0.5;
    static constexpr amrex::Real c20 =  c10;
    static constexpr amrex::Real c21 =  c11;
    static constexpr amrex::Real c22 =  c12;
    
    static constexpr amrex::Real s00 =  4./3.; // axx 
    static constexpr amrex::Real s01 = -4; // bxx
    static constexpr amrex::Real s02 =  8./3.; // cxx
    static constexpr amrex::Real s10 =  1; // ayy
    static constexpr amrex::Real s11 =  -2; // byy
    static constexpr amrex::Real s12 =  1; // cyy
    static constexpr amrex::Real s20 =  s10; // azz
    static constexpr amrex::Real s21 =  s11; // bzz
    static constexpr amrex::Real s22 =  s12; // czz

    static constexpr amrex::Real m00  =   0.25; // axy 
    static constexpr amrex::Real m01  =  -0.25; // bxy
    static constexpr amrex::Real m02  =   0.;   // cxy
    static constexpr amrex::Real m03  =   0.25; // dxy
    static constexpr amrex::Real m04  =  -0.25; // 3xy
    static constexpr amrex::Real m10  =  m00;   // axz 
    static constexpr amrex::Real m11  =  m01;   // bxz
    static constexpr amrex::Real m12  =  m02;   // cxz
    static constexpr amrex::Real m13  =  m03;   // dyz
    static constexpr amrex::Real m14  =  m04;   // eyz
    static constexpr amrex::Real m20  =  m00;   // ayz 
    static constexpr amrex::Real m21  =  m01;   // byz
    static constexpr amrex::Real m22  =  m02;   // cyz
    static constexpr amrex::Real m23  =  m03;   // dyz
    static constexpr amrex::Real m24  =  m04;   // eyz
};

struct StencilJLO
{
    static constexpr amrex::Real c00 =  0.5;
    static constexpr amrex::Real c01 =  0.0;
    static constexpr amrex::Real c02 = -0.5;
    static constexpr amrex::Real c10 =  1./3;
    static constexpr amrex::Real c11 =  1;
    static constexpr amrex::Real c12 = -4./3.;
    static constexpr amrex::Real c20 =  c00;
    static constexpr amrex::Real c21 =  c01;
    static constexpr amrex::Real c22 =  c02;

    static constexpr amrex::Real s00 =  1.0; // axx 
    static constexpr amrex::Real s01 = -2.0; // bxx
    static constexpr amrex::Real s02 =  1.0; // cxx
    static constexpr amrex::Real s10 =  4./3.; // ayy
    static constexpr amrex::Real s11 =  -4.; // byy
    static constexpr amrex::Real s12 =  8./3.; // cyy
    static constexpr amrex::Real s20 =  s00; // azz
    static constexpr amrex::Real s21 =  s01; // bzz
    static constexpr amrex::Real s22 =  s02; // czz

    static constexpr amrex::Real m00  =   0.25; // axy 
    static constexpr amrex::Real m01  =  -0.25; // bxy
    static constexpr amrex::Real m02  =   0.;   // cxy
    static constexpr amrex::Real m03  =   0.25; // dxy
    static constexpr amrex::Real m04  =  -0.25; // 3xy
    static constexpr amrex::Real m10  =  m00;   // axz 
    static constexpr amrex::Real m11  =  m01;   // bxz
    static constexpr amrex::Real m12  =  m02;   // cxz
    static constexpr amrex::Real m13  =  m03;   // dyz
    static constexpr amrex::Real m14  =  m04;   // eyz
    static constexpr amrex::Real m20  =  m00;   // ayz 
    static constexpr amrex::Real m21  =  m01;   // byz
    static constexpr amrex::Real m22  =  m02;   // cyz
    static constexpr amrex::Real m23  =  m03;   // dyz
    static constexpr amrex::Real m24  =  m04;   // eyz

};

struct StencilKLO
{
    static constexpr amrex::Real c00 =  0.5;
    static constexpr amrex::Real c01 =  0.0;
    static constexpr amrex::Real c02 = -0.5;
    static constexpr amrex::Real c10 =  c00;
    static constexpr amrex::Real c11 =  c01;
    static constexpr amrex::Real c12 =  c02;
    static constexpr amrex::Real c20 =  1./3.;
    static constexpr amrex::Real c21 =  1.;
    static constexpr amrex::Real c22 = -4./3.;
    
    static constexpr amrex::Real s00 =  1.0; // axx 
    static constexpr amrex::Real s01 = -2.0; // bxx
    static constexpr amrex::Real s02 =  1.0; // cxx
    static constexpr amrex::Real s10 =  s00; // ayy
    static constexpr amrex::Real s11 =  s01; // byy
    static constexpr amrex::Real s12 =  s02; // cyy
    static constexpr amrex::Real s20 =  4./3.; // azz
    static constexpr amrex::Real s21 =  -4.;   // bzz
    static constexpr amrex::Real s22 =  8./3.; // czz
     
    static constexpr amrex::Real m00  =   0.25; // axy 
    static constexpr amrex::Real m01  =  -0.25; // bxy
    static constexpr amrex::Real m02  =   0.;   // cxy
    static constexpr amrex::Real m03  =   0.25; // dxy
    static constexpr amrex::Real m04  =  -0.25; // 3xy
    static constexpr amrex::Real m10  =  m00;   // axz 
    static constexpr amrex::Real m11  =  m01;   // bxz
    static constexpr amrex::Real m12  =  m02;   // cxz
    static constexpr amrex::Real m13  =  m03;   // dyz
    static constexpr amrex::Real m14  =  m04;   // eyz
    static constexpr amrex::Real m20  =  m00;   // ayz 
    static constexpr amrex::Real m21  =  m01;   // byz
    static constexpr amrex::Real m22  =  m02;   // cyz
    static constexpr amrex::Real m23  =  m03;   // dyz
    static constexpr amrex::Real m24  =  m04;   // eyz

};

struct StencilIHI
{
    static constexpr amrex::Real c00 =  4./3.;
    static constexpr amrex::Real c01 =  -1;
    static constexpr amrex::Real c02 =  -1./3.;
    static constexpr amrex::Real c10 =  0.5;
    static constexpr amrex::Real c11 =  0.0;
    static constexpr amrex::Real c12 = -0.5;
    static constexpr amrex::Real c20 =  c10;
    static constexpr amrex::Real c21 =  c11;
    static constexpr amrex::Real c22 =  c12;

    static constexpr amrex::Real s00 = 8./3.; // axx 
    static constexpr amrex::Real s01 = -4.; // bxx
    static constexpr amrex::Real s02 = 4./3.; // cxx
    static constexpr amrex::Real s10 =  1.; // ayy
    static constexpr amrex::Real s11 = -2.; // byy
    static constexpr amrex::Real s12 =  1.; // cyy
    static constexpr amrex::Real s20 =  s10; // azz
    static constexpr amrex::Real s21 =  s11; // bzz
    static constexpr amrex::Real s22 =  s12; // czz


    static constexpr amrex::Real m00  =   0.25; // axy 
    static constexpr amrex::Real m01  =  -0.25; // bxy
    static constexpr amrex::Real m02  =   0.;   // cxy
    static constexpr amrex::Real m03  =   0.25; // dxy
    static constexpr amrex::Real m04  =  -0.25; // 3xy
    static constexpr amrex::Real m10  =  m00;   // axz 
    static constexpr amrex::Real m11  =  m01;   // bxz
    static constexpr amrex::Real m12  =  m02;   // cxz
    static constexpr amrex::Real m13  =  m03;   // dyz
    static constexpr amrex::Real m14  =  m04;   // eyz
    static constexpr amrex::Real m20  =  m00;   // ayz 
    static constexpr amrex::Real m21  =  m01;   // byz
    static constexpr amrex::Real m22  =  m02;   // cyz
    static constexpr amrex::Real m23  =  m03;   // dyz
    static constexpr amrex::Real m24  =  m04;   // eyz

};

struct StencilJHI
{
    static constexpr amrex::Real c00 =  0.5;
    static constexpr amrex::Real c01 =  0.0;
    static constexpr amrex::Real c02 = -0.5;
    static constexpr amrex::Real c10 =  4./3.;
    static constexpr amrex::Real c11 = -1.;
    static constexpr amrex::Real c12 = -1./3.;
    static constexpr amrex::Real c20 =  c00;
    static constexpr amrex::Real c21 =  c01;
    static constexpr amrex::Real c22 =  c02;
    
    static constexpr amrex::Real s00 =  1.0; // axx 
    static constexpr amrex::Real s01 = -2.0; // bxx
    static constexpr amrex::Real s02 =  1.0; // cxx
    static constexpr amrex::Real s10 =  8./3.; // ayy
    static constexpr amrex::Real s11 =  -4.; // byy
    static constexpr amrex::Real s12 =  4./3.; // cyy
    static constexpr amrex::Real s20 =  s00; // azz
    static constexpr amrex::Real s21 =  s01; // bzz
    static constexpr amrex::Real s22 =  s02; // czz

    static constexpr amrex::Real m00  =   0.25; // axy 
    static constexpr amrex::Real m01  =  -0.25; // bxy
    static constexpr amrex::Real m02  =   0.;   // cxy
    static constexpr amrex::Real m03  =   0.25; // dxy
    static constexpr amrex::Real m04  =  -0.25; // 3xy
    static constexpr amrex::Real m10  =  m00;   // axz 
    static constexpr amrex::Real m11  =  m01;   // bxz
    static constexpr amrex::Real m12  =  m02;   // cxz
    static constexpr amrex::Real m13  =  m03;   // dyz
    static constexpr amrex::Real m14  =  m04;   // eyz
    static constexpr amrex::Real m20  =  m00;   // ayz 
    static constexpr amrex::Real m21  =  m01;   // byz
    static constexpr amrex::Real m22  =  m02;   // cyz
    static constexpr amrex::Real m23  =  m03;   // dyz
    static constexpr amrex::Real m24  =  m04;   // eyz

};

struct StencilKHI
{
    static constexpr amrex::Real c00 =  0.5;
    static constexpr amrex::Real c01 =  0.0;
    static constexpr amrex::Real c02 = -0.5;
    static constexpr amrex::Real c10 =  c00;
    static constexpr amrex::Real c11 =  c01;
    static constexpr amrex::Real c12 =  c02;
    static constexpr amrex::Real c20 =  4./3.;
    static constexpr amrex::Real c21 = -1.;
    static constexpr amrex::Real c22 = -1./3.;

    static constexpr amrex::Real s00 =  1.0;   // axx 
    static constexpr amrex::Real s01 = -2.0;   // bxx
    static constexpr amrex::Real s02 =  1.0;   // cxx
    static constexpr amrex::Real s10 =  s00;   // ayy
    static constexpr amrex::Real s11 =  s01;   // byy
    static constexpr amrex::Real s12 =  s02;   // cyy
    static constexpr amrex::Real s20 =  8./3.; // azz
    static constexpr amrex::Real s21 =  -4.;   // bzz
    static constexpr amrex::Real s22 =  4./3.; // czz

    static constexpr amrex::Real m00  =   0.25; // axy 
    static constexpr amrex::Real m01  =  -0.25; // bxy
    static constexpr amrex::Real m02  =   0.;   // cxy
    static constexpr amrex::Real m03  =   0.25; // dxy
    static constexpr amrex::Real m04  =  -0.25; // 3xy
    static constexpr amrex::Real m10  =  m00;   // axz 
    static constexpr amrex::Real m11  =  m01;   // bxz
    static constexpr amrex::Real m12  =  m02;   // cxz
    static constexpr amrex::Real m13  =  m03;   // dyz
    static constexpr amrex::Real m14  =  m04;   // eyz
    static constexpr amrex::Real m20  =  m00;   // ayz 
    static constexpr amrex::Real m21  =  m01;   // byz
    static constexpr amrex::Real m22  =  m02;   // cyz
    static constexpr amrex::Real m23  =  m03;   // dyz
    static constexpr amrex::Real m24  =  m04;   // eyz

};



template <typename Stencil>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void gradient(
    int i,
    int j,
    int k,
    amrex::Real idx,
    amrex::Real idy,
    amrex::Real idz,
    amrex::Array4<amrex::Real const> const& phi,
    amrex::Array4<amrex::Real> const& gradphi,
    int ncomp) noexcept
{
    using namespace amrex;

    Real cp1, c, cm1;
    for (int icomp=0; icomp < ncomp; icomp++) {
        cp1 = Stencil::c00;
        c   = Stencil::c01;
        cm1 = Stencil::c02;
        gradphi(i,j,k,icomp*AMREX_SPACEDIM+0) = (cp1*phi(i+1,j,k,icomp) + c*phi(i,j,k,icomp) + cm1*phi(i-1,j,k,icomp)) * idx;
        cp1 = Stencil::c10;
        c   = Stencil::c11;
        cm1 = Stencil::c12;
        gradphi(i,j,k,icomp*AMREX_SPACEDIM+1) = (cp1*phi(i,j+1,k,icomp) + c*phi(i,j,k,icomp) + cm1*phi(i,j-1,k,icomp)) * idy;
        cp1 = Stencil::c20;
        c   = Stencil::c21;
        cm1 = Stencil::c22;
        gradphi(i,j,k,icomp*AMREX_SPACEDIM+2) = (cp1*phi(i,j,k+1,icomp) + c*phi(i,j,k,icomp) + cm1*phi(i,j,k-1,icomp)) * idz;
    }

}

template <typename Stencil>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void laplacian(
    int i,
    int j,
    int k,
    amrex::Real idx,
    amrex::Real idy,
    amrex::Real idz,
    amrex::Array4<amrex::Real const> const& phi,
    amrex::Array4<amrex::Real> const& Laplacephi,
    int ncomp) noexcept
{
    using namespace amrex;

    Real sp1, s, sm1, d2phidx2, d2phidy2,d2phidz2;
    for (int icomp=0; icomp < ncomp; icomp++) {
        sp1 = Stencil::s00;
        s   = Stencil::s01;
        sm1 = Stencil::s02;
        d2phidx2 = (sp1*phi(i+1,j,k,icomp) + s*phi(i,j,k,icomp) + sm1*phi(i-1,j,k,icomp))*idx*idx;
        sp1 = Stencil::s10;
        s   = Stencil::s11;
        sm1 = Stencil::s12;
        d2phidy2 = (sp1*phi(i,j+1,k,icomp) + s*phi(i,j,k,icomp) + sm1*phi(i,j-1,k,icomp))*idy*idy;
        sp1 = Stencil::s20;
        s   = Stencil::s21;
        sm1 = Stencil::s22;
        d2phidz2 = (sp1*phi(i,j,k+1,icomp) + s*phi(i,j,k,icomp) + sm1*phi(i,j,k-1,icomp))*idz*idz;
        
        Laplacephi(i,j,k,icomp) = d2phidx2+d2phidy2+d2phidz2; 
    }
}


template <typename Stencil>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real curvature(
    int i,
    int j,
    int k,
    amrex::Real idx,
    amrex::Real idy,
    amrex::Real idz,
    amrex::Array4<amrex::Real const> const& phi) noexcept
{
    using namespace amrex;

    Real cp1, c, cm1,sp1,s,sm1,mp1p1,mm1p1,m,mp1m1,mm1m1;

    cp1 = Stencil::c00;
    c   = Stencil::c01; 
    cm1 = Stencil::c02;
    const Real phix = (cp1*phi(i+1,j,k) + c*phi(i,j,k) + cm1*phi(i-1,j,k)) * idx;
    sp1 = Stencil::s00;
    s   = Stencil::s01; 
    sm1 = Stencil::s02;
    const Real phixx = (sp1*phi(i+1,j,k) + s*phi(i,j,k) + sm1*phi(i-1,j,k)) * idx*idx;
    
    cp1 = Stencil::c10;
    c   = Stencil::c11;
    cm1 = Stencil::c12;
    const Real phiy  = (cp1*phi(i,j+1,k) + c*phi(i,j,k) + cm1*phi(i,j-1,k)) * idy;
    sp1 = Stencil::s10;
    s   = Stencil::s11; 
    sm1 = Stencil::s12;
    const Real phiyy = (sp1*phi(i,j+1,k) + s*phi(i,j,k) + sm1*phi(i,j-1,k)) *idy*idy;
    mp1p1  = Stencil::m00; 
    mm1p1  = Stencil::m01; 
    m      = Stencil::m02;
    mp1m1  = Stencil::m03; 
    mm1m1  = Stencil::m04; 
    const Real phixy = (mp1p1*phi(i+1,j+1,k) + mm1p1*phi(i-1,j+1,k) + m*phi(i,j,k)+mp1m1*phi(i+1,j-1,k) + mm1m1*phi(i-1,j-1,k)) *idx*idy;
    
    cp1 = Stencil::c20;
    c   = Stencil::c21;
    cm1 = Stencil::c22;
    const Real phiz  = (cp1*phi(i,j,k+1) + c*phi(i,j,k) + cm1*phi(i,j,k-1)) * idz;
    sp1 = Stencil::s20;
    s   = Stencil::s21; 
    sm1 = Stencil::s22; 
    const Real phizz = (sp1*phi(i,j,k+1) + s*phi(i,j,k) + sm1*phi(i,j,k-1)) *idz*idz;
    mp1p1  = Stencil::m10; 
    mm1p1  = Stencil::m11; 
    m      = Stencil::m12;
    mp1m1  = Stencil::m13; 
    mm1m1  = Stencil::m14; 
    const Real phixz = (mp1p1*phi(i+1,j,k+1) + mm1p1*phi(i-1,j,k+1) + m*phi(i,j,k)+mp1m1*phi(i+1,j,k-1) + mm1m1*phi(i-1,j,k-1)) *idx*idz;
    mp1p1  = Stencil::m20; 
    mm1p1  = Stencil::m21; 
    m      = Stencil::m22;
    mp1m1  = Stencil::m23; 
    mm1m1  = Stencil::m24; 
    const Real phiyz = (mp1p1*phi(i,j+1,k+1) + mm1p1*phi(i,j-1,k+1) + m*phi(i,j,k)+mp1m1*phi(i,j+1,k-1) + mm1m1*phi(i,j-1,k-1)) *idy*idz; 
    
    Real curv= -(phix*phix*phiyy - 2.*phix*phiy*phixy + phiy*phiy*phixx + phix*phix*phizz
                -2.*phix*phiz*phixz + phiz*phiz*phixx + phiy*phiy*phizz - 2*phiy*phiz*phiyz 
                + phiz*phiz*phiyy)/std::pow(phix*phix+phiy*phiy+phiz*phiz,1.5);
       
    return curv;
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

template<typename FType>
void compute_gradient(FType& gradf, const Field& field);

template<typename FType>
void compute_laplacian(FType& laplacef, const Field& field);

template<typename FType>
void compute_curvature(FType& curvf, const Field& field);

} // namespace amr_wind
#endif
