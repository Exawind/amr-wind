#include <incflo_divop_conv.hpp>
#include <param_mod_F.H>

#include <cmath>
#include <limits>

namespace divop_conv_aux {

AMREX_GPU_HOST_DEVICE
Real 
interp_to_face_centroid_x(const int i,
                          const int j,
                          const int k,
                          Array4<Real> const& var,
                          const int n,
                          Array4<const Real> const& afrac,
                          Array4<const Real> const& cent,
                          const EBCellFlag& flag)
{
  Real result(0);

  Real frac_y(0), frac_z(0);

  if (afrac(i,j,k) == 0)
    result = 0;
  else if (afrac(i,j,k) == 1)
    result = var(i,j,k,n);
  else
  {
    if (cent(i,j,k,0) < 0)
    {
      frac_y = -1 * cent(i,j,k,0) * flag.isConnected({0,-1,0});
      if (cent(i,j,k,1) <= 0)
      {
        frac_z = -1 * cent(i,j,k,1) * flag.isConnected({0,0,-1});
        result = (1-frac_z) * (   frac_y * var(i,j-1,k  ,n)  +
                              (1-frac_y) * var(i,j  ,k  ,n)) +
                     frac_z * (   frac_y * var(i,j-1,k-1,n)  +
                              (1-frac_y) * var(i,j  ,k-1,n));
      }
      else
      {
        frac_z = cent(i,j,k,1) * flag.isConnected({0,0,1});
        result = (1-frac_z) * (   frac_y * var(i,j-1,k  ,n)  +
                              (1-frac_y) * var(i,j  ,k  ,n)) +
                     frac_z * (   frac_y * var(i,j-1,k+1,n)  +
                              (1-frac_y) * var(i,j  ,k+1,n));
      }
    }
    else
    {
      frac_y = cent(i,j,k,0) * flag.isConnected({0,1,0});
      if (cent(i,j,k,1) <= 0)
      {
        frac_z = -1 * cent(i,j,k,1) * flag.isConnected({0,0,-1});
        result = (1-frac_z) * (   frac_y * var(i,j+1,k  ,n)  +
                              (1-frac_y) * var(i,j  ,k  ,n)) +
                     frac_z * (   frac_y * var(i,j+1,k-1,n)  +
                              (1-frac_y) * var(i,j  ,k-1,n));
      }
      else
      {
        frac_z = cent(i,j,k,1) * flag.isConnected({0,0,1});
        result = (1-frac_z) * (   frac_y * var(i,j+1,k  ,n)  +
                              (1-frac_y) * var(i,j  ,k  ,n)) +
                     frac_z * (   frac_y * var(i,j+1,k+1,n)  +
                              (1-frac_y) * var(i,j  ,k+1,n));
      }
    }
  }

  return result;
}

AMREX_GPU_HOST_DEVICE
Real 
interp_to_face_centroid_y(const int i,
                          const int j,
                          const int k,
                          Array4<Real> const& var,
                          const int n,
                          Array4<const Real> const& afrac,
                          Array4<const Real> const& cent,
                          const EBCellFlag& flag)
{
  Real result(0);

  Real frac_x(0), frac_z(0);

  if (afrac(i,j,k) == 0)
    result = 0;
  else if (afrac(i,j,k) == 1)
    result = var(i,j,k,n);
  else
  {
    if (cent(i,j,k,0) < 0)
    {
      frac_x = -1 * cent(i,j,k,0) * flag.isConnected({-1,0,0});
      if (cent(i,j,k,1) <= 0)
      {
        frac_z = -1 * cent(i,j,k,1) * flag.isConnected({0,0,-1});
        result = (1-frac_z) * (   frac_x * var(i-1,j,k  ,n)  +
                              (1-frac_x) * var(i  ,j,k  ,n)) +
                     frac_z * (   frac_x * var(i-1,j,k-1,n)  +
                              (1-frac_x) * var(i  ,j,k-1,n));
      }
      else
      {
         frac_z =  cent(i,j,k,1) * flag.isConnected({0,0,1});
         result = (1-frac_z) * (   frac_x * var(i-1,j,k  ,n)  +
                               (1-frac_x) * var(i  ,j,k  ,n)) +
                      frac_z * (   frac_x * var(i-1,j,k+1,n)  +
                               (1-frac_x) * var(i  ,j,k+1,n));
      }
    }
    else
    {
      frac_x = cent(i,j,k,0) * flag.isConnected({1,0,0});
      if (cent(i,j,k,1) <= 0)
      {
        frac_z = -1 * cent(i,j,k,1) * flag.isConnected({0,0,-1});
        result = (1-frac_z) * (   frac_x * var(i+1,j,k  ,n)  +
                              (1-frac_x) * var(i  ,j,k  ,n)) +
                     frac_z * (   frac_x * var(i+1,j,k-1,n)  +
                              (1-frac_x) * var(i  ,j,k-1,n));
      }
      else
      {
        frac_z =  cent(i,j,k,1) * flag.isConnected({0,0,1});
        result = (1-frac_z) * (   frac_x * var(i+1,j,k  ,n)  +
                              (1-frac_x) * var(i  ,j,k  ,n)) +
                     frac_z * (   frac_x * var(i+1,j,k+1,n)  +
                              (1-frac_x) * var(i  ,j,k+1,n));
      }
    } 
  }

  return result;
}

AMREX_GPU_HOST_DEVICE
Real 
interp_to_face_centroid_z(const int i,
                          const int j,
                          const int k,
                          Array4<Real> const& var,
                          const int n,
                          Array4<const Real> const& afrac,
                          Array4<const Real> const& cent,
                          const EBCellFlag& flag)
{
  Real result(0);

  Real frac_x(0), frac_y(0);

  if (afrac(i,j,k) == 0)
    result = 0;
  else if (afrac(i,j,k) == 1)
    result = var(i,j,k,n);
  else
  {
    if (cent(i,j,k,0) < 0)
    {
      frac_x = -1 * cent(i,j,k,0) * flag.isConnected({-1,0,0});
      if (cent(i,j,k,1) <= 0)
      {
        frac_y = -1 * cent(i,j,k,1) * flag.isConnected({0,-1,0});
        result = (1-frac_y) * (   frac_x * var(i-1,j  ,k,n)  +
                              (1-frac_x) * var(i  ,j  ,k,n)) +
                     frac_y * (   frac_x * var(i-1,j-1,k,n)  +
                              (1-frac_x) * var(i  ,j-1,k,n));
      }
      else
      {
        frac_y =  cent(i,j,k,1) * flag.isConnected({0,1,0});
        result = (1-frac_y) * (   frac_x * var(i-1,j  ,k,n)  +
                              (1-frac_x) * var(i  ,j  ,k,n)) +
                     frac_y * (   frac_x * var(i-1,j+1,k,n)  +
                              (1-frac_x) * var(i  ,j+1,k,n));
      }
    }
    else
    {
      frac_x = cent(i,j,k,0) * flag.isConnected({1,0,0});
      if (cent(i,j,k,1) <= 0)
      {
        frac_y = -1 * cent(i,j,k,1) * flag.isConnected({0,-1,0});
        result = (1-frac_y) * (   frac_x * var(i+1,j  ,k,n)  +
                              (1-frac_x) * var(i  ,j  ,k,n)) +
                     frac_y * (   frac_x * var(i+1,j-1,k,n)  +
                              (1-frac_x) * var(i  ,j-1,k,n));
      }
      else
      {
        frac_y =  cent(i,j,k,1) * flag.isConnected({0,1,0});
        result = (1-frac_y) * (   frac_x * var(i+1,j  ,k,n)  +
                              (1-frac_x) * var(i  ,j  ,k,n)) +
                     frac_y * (   frac_x * var(i+1,j+1,k,n)  +
                              (1-frac_x) * var(i  ,j+1,k,n));
      }
    }
  }

  return result;
}

void
step2(const Box& grown1_bx,
      const Box& grown2_bx,
      MFIter* mfi,
      FArrayBox& optmp_fbx,
      FArrayBox& divc_fbx,
      FArrayBox& delm_fbx,
      const MultiFab* volfrac,
      FArrayBox& mask_fbx,
      const EBCellFlagFab& flags_fab)
{
  Array4<const EBCellFlag> const& flags = flags_fab.array();

  Array4<const Real> const& vfrac = volfrac->array(*mfi);

  Array4<Real> const& delm = delm_fbx.array();

  Array4<Real> const& optmp = optmp_fbx.array();
  Array4<Real> const& divc = divc_fbx.array();
  Array4<Real> const& mask = mask_fbx.array();

  AMREX_HOST_DEVICE_FOR_3D(grown2_bx, i, j, k,
  {
    optmp(i,j,k) = 0;
  });

  AMREX_HOST_DEVICE_FOR_3D(grown1_bx, i, j, k,
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
              vtot += epvfrac;
              divnc += epvfrac * divc(i+ii,j+jj,k+kk);
            }

      divnc /= vtot;
      epvfrac = vfrac(i,j,k);
      optmp(i,j,k) = (1 - vfrac(i,j,k)) * (divnc - divc(i,j,k));
      delm(i,j,k) = -1 * epvfrac * optmp(i,j,k);
    }
    else
      delm(i,j,k) = 0;
  });

  Gpu::synchronize();
}

void
step3(const Box& grown1_bx,
      MFIter* mfi,
      FArrayBox& optmp_fbx,
      FArrayBox& delm_fbx,
      const MultiFab* volfrac,
      FArrayBox& mask_fbx,
      const EBCellFlagFab& flags_fab)
{
  Array4<const EBCellFlag> const& flags = flags_fab.array();

  Array4<const Real> const& vfrac = volfrac->array(*mfi);

  Array4<Real> const& delm = delm_fbx.array();

  Array4<Real> const& optmp = optmp_fbx.array();
  Array4<Real> const& mask = mask_fbx.array();

  AMREX_HOST_DEVICE_FOR_3D(grown1_bx, i, j, k,
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
              Gpu::Atomic::Add(&optmp(i+ii,j+jj,k+kk),
                                delm(i,j,k) * wtot * mask(i,j,k));
#else
              optmp(i+ii,j+jj,k+kk) += 
                delm(i,j,k) * wtot * mask(i,j,k);
#endif
            }
    }
  });

  Gpu::synchronize();
}

} // end namespace divop_conv_aux

using namespace divop_conv_aux;

void
compute_divop_conv(
              Box& bx,
              MultiFab& conv,
              int conv_comp, int ncomp,
              MFIter* mfi,
              FArrayBox& fxfab,
              FArrayBox& fyfab,
              FArrayBox& fzfab,
              Array<const MultiCutFab*, AMREX_SPACEDIM>& areafrac,
              Array<const MultiCutFab*, AMREX_SPACEDIM>& facecent,
              const EBCellFlagFab& flags_fab,
              const MultiFab* volfrac,
              const MultiCutFab* bndrycent_fab,
              Box& domain,
              const int cyclic_x,
              const int cyclic_y,
              const int cyclic_z,
              const Real* dx)
{
  const amrex::Dim3 dom_low = amrex::lbound(domain);
  const amrex::Dim3 dom_high = amrex::ubound(domain);

  const Real i_dx (1/dx[0]), i_dy(1/dx[1]), i_dz(1/dx[2]);

  Array4<Real> const& divergence = conv.array(*mfi);

  Array4<Real> const& fx = fxfab.array();
  Array4<Real> const& fy = fyfab.array();
  Array4<Real> const& fz = fzfab.array();

  Array4<const Real> const& areafrac_x = areafrac[0]->array(*mfi);
  Array4<const Real> const& areafrac_y = areafrac[1]->array(*mfi);
  Array4<const Real> const& areafrac_z = areafrac[2]->array(*mfi);

  Array4<const Real> const& facecent_x = facecent[0]->array(*mfi);
  Array4<const Real> const& facecent_y = facecent[1]->array(*mfi);
  Array4<const Real> const& facecent_z = facecent[2]->array(*mfi);

  Array4<const EBCellFlag> const& flags = flags_fab.array();

  Array4<const Real> const& vfrac = volfrac->array(*mfi);

  const Real tolerance = std::numeric_limits<Real>::epsilon();

  if((std::abs(dx[0] - dx[1]) > tolerance) or 
     (std::abs(dx[0] - dx[2]) > tolerance) or
     (std::abs(dx[1] - dx[2]) > tolerance))
    amrex::Abort("Compute divop(): grid spacing must be uniform");

  const Box& grown1_bx = amrex::grow(bx,1);
  const Box& grown2_bx = amrex::grow(bx,2);

  FArrayBox delm_fbx(grown1_bx);

  FArrayBox optmp_fbx(grown2_bx);
  FArrayBox divc_fbx(grown2_bx);
  FArrayBox mask_fbx(grown2_bx);

  Array4<Real> const& optmp = optmp_fbx.array();
  Array4<Real> const& divc = divc_fbx.array();
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

  const Real my_huge = get_my_huge();

  //
  // We use the EB algorithm to compute the divergence at cell centers
  //
  for(unsigned int n(0); n < ncomp; ++n)
  {
    //
    // Step 1: compute conservative divergence on stencil (lo-2,hi-2)
    //

    AMREX_HOST_DEVICE_FOR_3D(grown2_bx, i, j, k,
    {
      if(flags(i,j,k).isCovered())
      {
        divc(i,j,k) = my_huge;
      }
      else if(flags(i,j,k).isSingleValued())
      {
        const EBCellFlag& flag = flags(i,j,k);

        Real fxp = interp_to_face_centroid_x(i+1, j, k, fx, n,
                                             areafrac_x, facecent_x, flag);
        Real fxm = interp_to_face_centroid_x(i  , j, k, fx, n,
                                             areafrac_x, facecent_x, flag);

        Real fyp = interp_to_face_centroid_y(i, j+1, k, fy, n,
                                             areafrac_y, facecent_y, flag);
        Real fym = interp_to_face_centroid_y(i, j  , k, fy, n,
                                             areafrac_y, facecent_y, flag);
        
        Real fzp = interp_to_face_centroid_z(i, j, k+1, fz, n,
                                             areafrac_z, facecent_z, flag);
        Real fzm = interp_to_face_centroid_z(i, j, k  , fz, n,
                                             areafrac_z, facecent_z, flag);

        divc(i,j,k) = ((fxp*areafrac_x(i+1,j,k) - (fxm*areafrac_x(i,j,k))) * i_dx +
                       (fyp*areafrac_y(i,j+1,k) - (fym*areafrac_y(i,j,k))) * i_dy +
                       (fzp*areafrac_z(i,j,k+1) - (fzm*areafrac_z(i,j,k))) * i_dz) / vfrac(i,j,k);
      }
      else 
      {
        divc(i,j,k) = ((fx(i+1,j,k,n) - fx(i,j,k,n)) * i_dx +
                       (fy(i,j+1,k,n) - fy(i,j,k,n)) * i_dy +
                       (fz(i,j,k+1,n) - fz(i,j,k,n)) * i_dz);
       
      }
    });

    //
    // Step 2: compute delta M (mass gain or loss) on (lo-1,lo+1)
    //
    step2(grown1_bx, grown2_bx, mfi, optmp_fbx, divc_fbx, delm_fbx, 
        volfrac, mask_fbx, flags_fab);

    //
    // Step 3: redistribute excess/loss of mass
    //
    step3(grown1_bx, mfi, optmp_fbx, delm_fbx, 
        volfrac, mask_fbx, flags_fab);

    //
    // Resume the correct sign, AKA return the negative
    //
    AMREX_HOST_DEVICE_FOR_3D(bx, i, j, k,
    {
      divergence(i,j,k,conv_comp+n) = divc(i,j,k) + optmp(i,j,k);
    });

  }

  Gpu::synchronize();

}
