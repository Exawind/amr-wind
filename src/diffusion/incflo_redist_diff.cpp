#include <incflo_redist_diff.hpp>
#include <param_mod_F.H>

#include <cmath>
#include <limits>

namespace redist_diff_aux {

void
step2(const Box& grown1_bx,
      const Box& grown2_bx,
      const int n,
      MFIter* mfi,
      MultiFab& divtau_aux,
      FArrayBox& delm_fbx,
      FArrayBox& optmp_fbx,
      FArrayBox& mask_fbx,
      const MultiFab* volfrac,
      const EBCellFlagFab& flags_fab)
{
  Array4<Real> const& divc = divtau_aux.array(*mfi);

  Array4<Real> const& delm = delm_fbx.array();
  Array4<Real> const& optmp = optmp_fbx.array();
  Array4<Real> const& mask = mask_fbx.array();

  Array4<const Real> const& vfrac = volfrac->array(*mfi);

  Array4<const EBCellFlag> const& flags = flags_fab.array();

  AMREX_FOR_3D(grown1_bx, i, j, k,
  {
    if(flags(i,j,k).isSingleValued())
    {
      Real divnc = 0;
      Real vtot = 0;

      Real my_vfrac = 0;

      for(int ii(-1); ii <= 1; ii++)
        for(int jj(-1); jj <= 1; jj++)
          for(int kk(-1); kk <= 1; kk++)
            if((ii != 0 or jj != 0 or kk != 0) and 
                (flags(i,j,k).isConnected(ii,jj,kk) == 1))
            {
              my_vfrac = vfrac(i+ii,j+jj,k+kk) * mask(i+ii,j+jj,k+kk);
              vtot += my_vfrac;
              divnc += my_vfrac * divc(i+ii,j+jj,k+kk,n);
            }

      divnc /= vtot;
      my_vfrac = vfrac(i,j,k);
      optmp(i,j,k) = (1 - vfrac(i,j,k)) * (divnc - divc(i,j,k,n));
      delm(i,j,k) = (-1 * my_vfrac * optmp(i,j,k)) * mask(i,j,k);
    }
    else
      delm(i,j,k) = 0;
  });

  Gpu::synchronize();
}

void
step3(const Box& grown1_bx,
      MFIter* mfi,
      FArrayBox& delm_fbx,
      FArrayBox& optmp_fbx,
      FArrayBox& mask_fbx,
      const MultiFab* volfrac,
      const EBCellFlagFab& flags_fab)
{

  Array4<Real> const& delm = delm_fbx.array();
  Array4<Real> const& optmp = optmp_fbx.array();
  Array4<Real> const& mask = mask_fbx.array();

  Array4<const Real> const& vfrac = volfrac->array(*mfi);

  Array4<const EBCellFlag> const& flags = flags_fab.array();

  AMREX_FOR_3D(grown1_bx, i, j, k,
  {
    if(flags(i,j,k).isSingleValued())
    {
      Real wtot = 0;
      
      for(int ii(-1); ii <= 1; ii++)
        for(int jj(-1); jj <= 1; jj++)
          for(int kk(-1); kk <= 1; kk++)
            if((ii != 0 or jj != 0 or kk != 0) and
                (flags(i,j,k).isConnected(ii,jj,kk) == 1))
            {
              wtot += vfrac(i+ii,j+jj,k+kk) * mask(i+ii,j+jj,k+kk);
            }

      wtot = 1/wtot;
      
      for(int ii(-1); ii <= 1; ii++)
        for(int jj(-1); jj <= 1; jj++)
          for(int kk(-1); kk <= 1; kk++)
            if((ii != 0 or jj != 0 or kk != 0) and
                (flags(i,j,k).isConnected(ii,jj,kk) == 1))
            {
              // Note: delm has already been multiplied by mask
#ifdef AMREX_USE_CUDA
              Gpu::Atomic::Add(&optmp(i+ii,j+jj,k+kk), delm(i,j,k) * wtot);
#else
              optmp(i+ii,j+jj,k+kk) += delm(i,j,k) * wtot;
#endif
            }
    }
  });

  Gpu::synchronize();
}

} // end namespace redist_diff_aux

using namespace redist_diff_aux;

void
compute_redist_diff(Box& bx,
                    MultiFab& divtau,
                    MultiFab& divtau_aux,
                    MFIter* mfi,
                    const EBCellFlagFab& flags_fab,
                    const MultiFab* volfrac,
                    const int cyclic_x,
                    const int cyclic_y,
                    const int cyclic_z,
                    Box& domain,
                    const Real* dx)
{
  const amrex::Dim3 dom_low = amrex::lbound(domain);
  const amrex::Dim3 dom_high = amrex::ubound(domain);

  Array4<Real> const& divergence = divtau.array(*mfi);
  Array4<Real> const& divc = divtau_aux.array(*mfi);

  const Box& grown1_bx = amrex::grow(bx,1);
  const Box& grown2_bx = amrex::grow(bx,2);

  FArrayBox delm_fbx(grown1_bx);
  FArrayBox optmp_fbx(grown2_bx);
  FArrayBox mask_fbx(grown2_bx);

  Array4<Real> const& optmp = optmp_fbx.array();
  Array4<Real> const& mask = mask_fbx.array();

  //
  // Array "mask" is used to sever the link to ghost cells when the BCs are not
  // periodic
  // It is set to 1 when a cell can be used in computations, 0 otherwise
  //
  AMREX_FOR_3D(grown2_bx, i, j, k,
  {
    if(((not cyclic_x) and (i < dom_low.x or i > dom_high.x)) or
       ((not cyclic_y) and (j < dom_low.y or j > dom_high.y)) or
       ((not cyclic_z) and (k < dom_low.z or k > dom_high.z)))
      mask(i,j,k) = 0;
    else
      mask(i,j,k) = 1;
  });

  //
  // Here we do the redistribution steps
  //
  for(unsigned int n(0); n < 3; ++n)
  {
    // Set this to zero here
    optmp_fbx.setVal(0.0);
 
    //
    // Step 2: compute delta M (mass gain or loss) on (lo-1,lo+1)
    //
    step2(grown1_bx, grown2_bx, n, mfi, divtau_aux, delm_fbx, optmp_fbx,
        mask_fbx, volfrac, flags_fab);

    //
    // Step 3: redistribute excess/loss of mass
    //
    step3(grown1_bx, mfi, delm_fbx, optmp_fbx, mask_fbx, volfrac, flags_fab);

    //
    // Resume the correct sign, AKA return the negative
    //
    AMREX_FOR_3D(bx, i, j, k,
    {
      divergence(i,j,k,n) = divc(i,j,k,n) + optmp(i,j,k);
    });
  }

  Gpu::synchronize();
}
