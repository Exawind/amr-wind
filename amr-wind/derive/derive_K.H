#ifndef DERIVE_K_H_
#define DERIVE_K_H_

#include "amr-wind/core/Field.H"
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
    
    static constexpr amrex::Real s00 =  4./3.;  // axx 
    static constexpr amrex::Real s01 = -4;      // bxx
    static constexpr amrex::Real s02 =  8./3.;  // cxx
    static constexpr amrex::Real s10 =  1;      // ayy
    static constexpr amrex::Real s11 =  -2;     // byy
    static constexpr amrex::Real s12 =  1;      // cyy
    static constexpr amrex::Real s20 =  s10;    // azz
    static constexpr amrex::Real s21 =  s11;    // bzz
    static constexpr amrex::Real s22 =  s12;    // czz

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

    static constexpr amrex::Real s00 =  1.0;    // axx 
    static constexpr amrex::Real s01 = -2.0;    // bxx
    static constexpr amrex::Real s02 =  1.0;    // cxx
    static constexpr amrex::Real s10 =  4./3.;  // ayy
    static constexpr amrex::Real s11 =  -4.;    // byy
    static constexpr amrex::Real s12 =  8./3.;  // cyy
    static constexpr amrex::Real s20 =  s00;    // azz
    static constexpr amrex::Real s21 =  s01;    // bzz
    static constexpr amrex::Real s22 =  s02;    // czz

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
    
    static constexpr amrex::Real s00 =  1.0;    // axx 
    static constexpr amrex::Real s01 = -2.0;    // bxx
    static constexpr amrex::Real s02 =  1.0;    // cxx
    static constexpr amrex::Real s10 =  s00;    // ayy
    static constexpr amrex::Real s11 =  s01;    // byy
    static constexpr amrex::Real s12 =  s02;    // cyy
    static constexpr amrex::Real s20 =  4./3.;  // azz
    static constexpr amrex::Real s21 =  -4.;    // bzz
    static constexpr amrex::Real s22 =  8./3.;  // czz
     
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

    static constexpr amrex::Real s00 = 8./3.;   // axx 
    static constexpr amrex::Real s01 = -4.;     // bxx
    static constexpr amrex::Real s02 = 4./3.;   // cxx
    static constexpr amrex::Real s10 =  1.;     // ayy
    static constexpr amrex::Real s11 = -2.;     // byy
    static constexpr amrex::Real s12 =  1.;     // cyy
    static constexpr amrex::Real s20 =  s10;    // azz
    static constexpr amrex::Real s21 =  s11;    // bzz
    static constexpr amrex::Real s22 =  s12;    // czz

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
    
    static constexpr amrex::Real s00 =  1.0;    // axx 
    static constexpr amrex::Real s01 = -2.0;    // bxx
    static constexpr amrex::Real s02 =  1.0;    // cxx
    static constexpr amrex::Real s10 =  8./3.;  // ayy
    static constexpr amrex::Real s11 =  -4.;    // byy
    static constexpr amrex::Real s12 =  4./3.;  // cyy
    static constexpr amrex::Real s20 =  s00;    // azz
    static constexpr amrex::Real s21 =  s01;    // bzz
    static constexpr amrex::Real s22 =  s02;    // czz
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
    amrex::Array4<amrex::Real const> const& gradphi) noexcept
{
    using namespace amrex;

    Real cp1, c, cm1;

    cp1 = Stencil::c00;
    c   = Stencil::c01; 
    cm1 = Stencil::c02;
    const Real phix  = gradphi(i,j,k,0);
    const Real phixx = (cp1*gradphi(i+1,j,k,0) + c*gradphi(i,j,k,0) + cm1*gradphi(i-1,j,k,0)) * idx;
    
    cp1 = Stencil::c10;
    c   = Stencil::c11;
    cm1 = Stencil::c12;
    const Real phiy  = gradphi(i,j,k,1);
    const Real phixy = (cp1*gradphi(i,j+1,k,0) + c*gradphi(i,j,k,0) + cm1*gradphi(i,j-1,k,0))*idy;
    const Real phiyy = (cp1*gradphi(i,j+1,k,1) + c*gradphi(i,j,k,1) + cm1*gradphi(i,j-1,k,1))*idy;
    
    cp1 = Stencil::c20;
    c   = Stencil::c21;
    cm1 = Stencil::c22;
    const Real phiz  = gradphi(i,j,k,2);
    const Real phixz = (cp1*gradphi(i,j,k+1,0) + c*gradphi(i,j,k,0) + cm1*gradphi(i,j,k-1,0)) *idz;
    const Real phiyz = (cp1*gradphi(i,j,k+1,1) + c*gradphi(i,j,k,1) + cm1*gradphi(i,j,k-1,1)) *idz;
    const Real phizz = (cp1*gradphi(i,j,k+1,2) + c*gradphi(i,j,k,2) + cm1*gradphi(i,j,k-1,2)) *idz;
    
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

/** Compute strain-rate for a given field
 *  \ingroup field_ops
 *
 *  \tparam FType Field or ScratchField
 *  \param field [out] The computed magnitude of strain-rate
 *  \param field [in] The velocity field
 */
template<typename FType>
void compute_strainrate(FType& field, const Field& velocity);

/** Compute gradient of a scalar or vector field
 *  \ingroup field_ops
 *
 *  \f$\nabla \mathbf{f}\f$
 *
 *  \tparam FType Field or ScratchField
 *  \param gradf [out] The gradient of the input field
 *  \param field [in] Field on which gradient is to be computed
 */
template<typename FType>
void compute_gradient(FType& gradf, const Field& field);

/** Compute laplacian of a scalar or vector field
 *  \ingroup field_ops
 *
 *  \f$\nabla^2 \mathbf{f}\f$
 *
 *  \tparam FType Field or ScratchField
 *  \param gradf [out] The gradient of the input field
 *  \param field [in] Field on which gradient is to be computed
 */
template<typename FType>
void compute_laplacian(FType& laplacef, const Field& field);

template<typename FType>
void compute_curvature(FType& curvf, const Field& field);

/** Normalize a given vector field
 *  \ingroup field_ops
 *
 *  \f$ \mathbf{f} = \mathbf{f} / |\mathbf{f}|\f$
 *  \tparam FType Field or ScratchField
 *  \param field [inout] Field to be normalized
 */
template<typename FType>
void normalize_field(FType& field);

} // namespace amr_wind
#endif
