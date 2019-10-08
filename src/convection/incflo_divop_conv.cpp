#include <incflo_divop_conv.hpp>
#include <param_mod_F.H>

#include <cmath>
#include <limits>

namespace divop_conv_aux {

void
step2(const Box& grown1_bx,
      const Box& grown2_bx,
      MFIter* mfi,
      FArrayBox& optmp_fbx,
      FArrayBox& divc_fbx,
      FArrayBox& delm_fbx,
      const MultiFab* volfrac,
      FArrayBox& mask_fbx,
      const EBCellFlagFab& flags_fab,
      const int icomp, const int ncomp)
{
  Array4<const EBCellFlag> const& flags = flags_fab.array();

  Array4<const Real> const& vfrac = volfrac->array(*mfi);

  Array4<Real> const& delm = delm_fbx.array();

  Array4<Real> const& optmp = optmp_fbx.array();
  Array4<Real> const& divc = divc_fbx.array();
  Array4<Real> const& mask = mask_fbx.array();

  AMREX_HOST_DEVICE_FOR_4D(grown2_bx, ncomp, i, j, k, n,
  {
    optmp(i,j,k,n) = 0;
  });

  AMREX_HOST_DEVICE_FOR_4D(grown1_bx, ncomp, i, j, k, n,
  {
    if(flags(i,j,k).isSingleValued())
    {
      Real divnc = 0;
      Real vtot = 0;
      Real epvfrac = 0;

      for(int ii(-1); ii <= 1; ii++)
        for(int jj(-1); jj <= 1; jj++)
          for(int kk(-1); kk <= 1; kk++)
            if((ii != 0 or jj != 0 or kk != 0) and 
                (flags(i,j,k).isConnected({ii,jj,kk}) == 1))
            {
              epvfrac = vfrac(i+ii,j+jj,k+kk) * mask(i+ii,j+jj,k+kk);
              vtot  += epvfrac;
              divnc += epvfrac * divc(i+ii,j+jj,k+kk,n);
            }

      divnc /= vtot;

      optmp(i,j,k,n) =  (1 - vfrac(i,j,k)) * (divnc - divc(i,j,k,n));
       delm(i,j,k,n) = -(    vfrac(i,j,k)) * optmp(i,j,k,n);
    }
    else
      delm(i,j,k,n) = 0;
  });

  AMREX_HOST_DEVICE_FOR_4D(grown1_bx, ncomp, i, j, k, n,
  {
    if(flags(i,j,k).isSingleValued())
    {
      Real wtot = 0;
      
      for(int ii(-1); ii <= 1; ii++)
        for(int jj(-1); jj <= 1; jj++)
          for(int kk(-1); kk <= 1; kk++)
            if((ii != 0 or jj != 0 or kk != 0) and
                (flags(i,j,k).isConnected({ii,jj,kk}) == 1))
            {
              wtot += vfrac(i+ii,j+jj,k+kk) * mask(i+ii,j+jj,k+kk);
            }

      wtot = 1/wtot;

      
      for(int ii(-1); ii <= 1; ii++)
        for(int jj(-1); jj <= 1; jj++)
          for(int kk(-1); kk <= 1; kk++)
            if((ii != 0 or jj != 0 or kk != 0) and
                (flags(i,j,k).isConnected({ii,jj,kk}) == 1))
            {
#ifdef AMREX_USE_CUDA
              Gpu::Atomic::Add(&optmp(i+ii,j+jj,k+kk,n), delm(i,j,k,n) * wtot * mask(i+ii,j+jj,k+kk);
#else
              optmp(i+ii,j+jj,k+kk,n) += delm(i,j,k,n) * wtot * mask(i+ii,j+jj,k+kk);
#endif
            }
    }
  });

  Gpu::streamSynchronize();
}

} // end namespace divop_conv_aux

using namespace divop_conv_aux;

void
incflo_apply_eb_redistribution ( Box& bx,
                               MultiFab& conv,
                               MultiFab& divc,
                               MFIter* mfi,
                               const int icomp,
                               const int ncomp,
                               const EBCellFlagFab& flags_fab,
                               const MultiFab* volfrac,
                               Box& domain,
                               const int cyclic_x,
                               const int cyclic_y,
                               const int cyclic_z,
                               const Real* dx)
{
  const amrex::Dim3 dom_low = amrex::lbound(domain);
  const amrex::Dim3 dom_high = amrex::ubound(domain);

  Array4<Real> const& divergence = conv.array(*mfi);

  const Real tolerance = std::numeric_limits<Real>::epsilon();

  if((std::abs(dx[0] - dx[1]) > tolerance) or 
     (std::abs(dx[0] - dx[2]) > tolerance) or
     (std::abs(dx[1] - dx[2]) > tolerance))
    amrex::Abort("Compute divop(): grid spacing must be uniform");

  const Box& grown1_bx = amrex::grow(bx,1);
  const Box& grown2_bx = amrex::grow(bx,2);

  FArrayBox  delm_fbx(grown1_bx,ncomp);
  FArrayBox  optmp_fbx(grown2_bx,ncomp);
  FArrayBox  mask_fbx(grown2_bx); 

  FArrayBox& divc_fbx =  divc[*mfi];    

  Array4<Real> const& optmp = optmp_fbx.array();
  Array4<Real> const& mask = mask_fbx.array();

  //
  // Array "mask" is used to sever the link to ghost cells when the BCs are not
  // periodic
  // It is set to 1 when a cell can be used in computations, 0 otherwise
  //
  AMREX_HOST_DEVICE_FOR_3D(grown2_bx, i, j, k,
  {
    if(((not cyclic_x) and (i < dom_low.x or i > dom_high.x)) or
       ((not cyclic_y) and (j < dom_low.y or j > dom_high.y)) or
       ((not cyclic_z) and (k < dom_low.z or k > dom_high.z)))
      mask(i,j,k) = 0;
    else
      mask(i,j,k) = 1;
  });

  //
  // Step 2: compute delta M (mass gain or loss) on (lo-1,lo+1)
  //
  step2(grown1_bx, grown2_bx, mfi, optmp_fbx, divc_fbx, delm_fbx, 
        volfrac, mask_fbx, flags_fab, icomp, ncomp);

  //
  // Resume the correct sign, AKA return the negative
  //
  Array4<Real> const& divcarr = divc.array(*mfi);

  AMREX_HOST_DEVICE_FOR_4D(bx, ncomp, i, j, k, n, 
  {
     divergence(i,j,k,icomp+n) = divcarr(i,j,k,n) + optmp(i,j,k,n);
  });

  Gpu::streamSynchronize();
}
