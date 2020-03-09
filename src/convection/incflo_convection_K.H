#ifndef INCFLO_CONVECTION_K_H_
#define INCFLO_CONVECTION_K_H_

#ifdef AMREX_USE_EB
#include <AMReX_EBFArrayBox.H>
#include <AMReX_EBCellFlag.H>
#else
#include <AMReX_FArrayBox.H>
#endif

namespace {

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real incflo_xslope (int i, int j, int k, int n,
                           amrex::Array4<amrex::Real const> const& vcc) noexcept
{
    amrex::Real dl = 2.0*(vcc(i  ,j,k,n) - vcc(i-1,j,k,n));
    amrex::Real dr = 2.0*(vcc(i+1,j,k,n) - vcc(i  ,j,k,n));
    amrex::Real dc = 0.5*(vcc(i+1,j,k,n) - vcc(i-1,j,k,n));
    amrex::Real slope = amrex::min(std::abs(dl),std::abs(dc),std::abs(dr));
    slope = (dr*dl > 0.0) ? slope : 0.0;
    return (dc > 0.0) ? slope : -slope;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real incflo_xslope_extdir (int i, int j, int k, int n,
                                  amrex::Array4<amrex::Real const> const& vcc,
                                  bool edlo, bool edhi, int domlo, int domhi) noexcept
{
    amrex::Real dl = 2.0*(vcc(i  ,j,k,n) - vcc(i-1,j,k,n));
    amrex::Real dr = 2.0*(vcc(i+1,j,k,n) - vcc(i  ,j,k,n));
    amrex::Real dc = 0.5*(vcc(i+1,j,k,n) - vcc(i-1,j,k,n));
    if (edlo and i == domlo) {
        dc = (vcc(i+1,j,k,n)+3.0*vcc(i,j,k,n)-4.0*vcc(i-1,j,k,n))/3.0;
    } else if (edhi and i == domhi) {
        dc = (4.0*vcc(i+1,j,k,n)-3.0*vcc(i,j,k,n)-vcc(i-1,j,k,n))/3.0;
    }
    amrex::Real slope = amrex::min(std::abs(dl),std::abs(dc),std::abs(dr));
    slope = (dr*dl > 0.0) ? slope : 0.0;
    return (dc > 0.0) ? slope : -slope;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real incflo_yslope (int i, int j, int k, int n,
                           amrex::Array4<amrex::Real const> const& vcc) noexcept
{
    amrex::Real dl = 2.0*(vcc(i,j  ,k,n) - vcc(i,j-1,k,n));
    amrex::Real dr = 2.0*(vcc(i,j+1,k,n) - vcc(i,j  ,k,n));
    amrex::Real dc = 0.5*(vcc(i,j+1,k,n) - vcc(i,j-1,k,n));
    amrex::Real slope = amrex::min(std::abs(dl),std::abs(dc),std::abs(dr));
    slope = (dr*dl > 0.0) ? slope : 0.0;
    return (dc > 0.0) ? slope : -slope;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real incflo_yslope_extdir (int i, int j, int k, int n,
                                  amrex::Array4<amrex::Real const> const& vcc,
                                  bool edlo, bool edhi, int domlo, int domhi) noexcept
{
    amrex::Real dl = 2.0*(vcc(i,j  ,k,n) - vcc(i,j-1,k,n));
    amrex::Real dr = 2.0*(vcc(i,j+1,k,n) - vcc(i,j  ,k,n));
    amrex::Real dc = 0.5*(vcc(i,j+1,k,n) - vcc(i,j-1,k,n));
    if (edlo and j == domlo) {
        dc = (vcc(i,j+1,k,n)+3.0*vcc(i,j,k,n)-4.0*vcc(i,j-1,k,n))/3.0;
    } else if (edhi and j == domhi) {
        dc = (4.0*vcc(i,j+1,k,n)-3.0*vcc(i,j,k,n)-vcc(i,j-1,k,n))/3.0;
    }
    amrex::Real slope = amrex::min(std::abs(dl),std::abs(dc),std::abs(dr));
    slope = (dr*dl > 0.0) ? slope : 0.0;
    return (dc > 0.0) ? slope : -slope;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real incflo_zslope (int i, int j, int k, int n,
                           amrex::Array4<amrex::Real const> const& vcc) noexcept
{
    amrex::Real dl = 2.0*(vcc(i,j,k  ,n) - vcc(i,j,k-1,n));
    amrex::Real dr = 2.0*(vcc(i,j,k+1,n) - vcc(i,j,k  ,n));
    amrex::Real dc = 0.5*(vcc(i,j,k+1,n) - vcc(i,j,k-1,n));
    amrex::Real slope = amrex::min(std::abs(dl),std::abs(dc),std::abs(dr));
    slope = (dr*dl > 0.0) ? slope : 0.0;
    return (dc > 0.0) ? slope : -slope;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real incflo_zslope_extdir (int i, int j, int k, int n,
                                  amrex::Array4<amrex::Real const> const& vcc,
                                  bool edlo, bool edhi, int domlo, int domhi) noexcept
{
    amrex::Real dl = 2.0*(vcc(i,j,k  ,n) - vcc(i,j,k-1,n));
    amrex::Real dr = 2.0*(vcc(i,j,k+1,n) - vcc(i,j,k  ,n));
    amrex::Real dc = 0.5*(vcc(i,j,k+1,n) - vcc(i,j,k-1,n));
    if (edlo and k == domlo) {
        dc = (vcc(i,j,k+1,n)+3.0*vcc(i,j,k,n)-4.0*vcc(i,j,k-1,n))/3.0;
    } else if (edhi and k == domhi) {
        dc = (4.0*vcc(i,j,k+1,n)-3.0*vcc(i,j,k,n)-vcc(i,j,k-1,n))/3.0;
    }
    amrex::Real slope = amrex::min(std::abs(dl),std::abs(dc),std::abs(dr));
    slope = (dr*dl > 0.0) ? slope : 0.0;
    return (dc > 0.0) ? slope : -slope;
}

#ifdef AMREX_USE_EB

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::GpuArray<amrex::Real,3>
incflo_slopes_eb (int i, int j, int k, int n,
                  amrex::Array4<amrex::Real const> const& state,
                  amrex::Array4<amrex::Real const> const& ccent,
                  amrex::Array4<amrex::EBCellFlag const> const& flag) noexcept
{
    amrex::Real A[27][3];
    amrex::Real du[27];

    {
      int lc=0;
      for(int kk(-1); kk<=1; kk++){
        for(int jj(-1); jj<=1; jj++){
          for(int ii(-1); ii<=1; ii++){

            if( flag(i,j,k).isConnected(ii,jj,kk) and
                not (ii==0 and jj==0 and kk==0)) {

            // Not multplying by dx to be consistent with how the
            // slope is stored. Also not including the global shift
            // wrt plo or i,j,k. We only need relative distance.

              A[lc][0] = ii + ccent(i+ii,j+jj,k+kk,0) - ccent(i,j,k,0);
              A[lc][1] = jj + ccent(i+ii,j+jj,k+kk,1) - ccent(i,j,k,1);
              A[lc][2] = kk + ccent(i+ii,j+jj,k+kk,2) - ccent(i,j,k,2);

              du[lc] = state(i+ii,j+jj,k+kk,n) - state(i,j,k,n);

            } else {

              A[lc][0] = 0.0;
              A[lc][1] = 0.0;
              A[lc][2] = 0.0;

              du[lc] = 0.0;
            }

            lc++;
          }
        }
      }
    }

    amrex::Real AtA[3][3];
    amrex::Real Atb[3];

    for(int jj(0); jj<3; ++jj){
      for(int ii(0); ii<3; ++ii){
        AtA[ii][jj] = 0.0;
      }
      Atb[jj] = 0.0;
    }

    {

      for(int lc(0); lc<27; ++lc){
        AtA[0][0] += A[lc][0]* A[lc][0];
        AtA[0][1] += A[lc][0]* A[lc][1];
        AtA[0][2] += A[lc][0]* A[lc][2];
        AtA[1][1] += A[lc][1]* A[lc][1];
        AtA[1][2] += A[lc][1]* A[lc][2];
        AtA[2][2] += A[lc][2]* A[lc][2];

        Atb[0] += A[lc][0]*du[lc];
        Atb[1] += A[lc][1]*du[lc];
        Atb[2] += A[lc][2]*du[lc];
      }
    }

    // Fill in symmetric
    AtA[1][0] = AtA[0][1];
    AtA[2][0] = AtA[0][2];
    AtA[2][1] = AtA[1][2];

    amrex::Real detAtA =
      AtA[0][0]*(AtA[1][1]*AtA[2][2] - AtA[1][2]*AtA[1][2]) -
      AtA[0][1]*(AtA[1][0]*AtA[2][2] - AtA[1][2]*AtA[2][0]) +
      AtA[0][2]*(AtA[1][0]*AtA[2][1] - AtA[1][1]*AtA[2][0]);

    amrex::Real xs = 0.0;
    amrex::Real ys = 0.0;
    amrex::Real zs = 0.0;

    // X direction
    if(flag(i  ,j,k).isSingleValued() or
      (flag(i-1,j,k).isSingleValued() or not flag(i,j,k).isConnected(-1,0,0)) or
      (flag(i+1,j,k).isSingleValued() or not flag(i,j,k).isConnected( 1,0,0))) 
    {

      amrex::Real detAtA_x =
        Atb[0]   *(AtA[1][1]*AtA[2][2] - AtA[1][2]*AtA[1][2]) -
        AtA[0][1]*(Atb[1] *  AtA[2][2] - AtA[1][2]*Atb[2]   ) +
        AtA[0][2]*(Atb[1] *  AtA[2][1] - AtA[1][1]*Atb[2]   );

      // Slope at centroid of (i,j,k)
      xs = detAtA_x / detAtA;

    } else {

      amrex::Real du_xl = 2.0*(state(i  ,j,k,n) - state(i-1,j,k,n));
      amrex::Real du_xr = 2.0*(state(i+1,j,k,n) - state(i  ,j,k,n));
      amrex::Real du_xc = 0.5*(state(i+1,j,k,n) - state(i-1,j,k,n));

      amrex::Real xslope = amrex::min(std::abs(du_xl), std::abs(du_xc), std::abs(du_xr));
      xslope = (du_xr*du_xl > 0.0) ? xslope : 0.0;
      xs = (du_xc > 0.0) ? xslope : -xslope;

    }

    // Y direction
    if(flag(i,j  ,k).isSingleValued() or
      (flag(i,j-1,k).isSingleValued() or not flag(i,j,k).isConnected(0,-1,0)) or
      (flag(i,j+1,k).isSingleValued() or not flag(i,j,k).isConnected(0, 1,0))) 
    {

      amrex::Real detAtA_y =
        AtA[0][0]*(Atb[1]  * AtA[2][2] - AtA[1][2]*Atb[2]   ) -
        Atb[0] *  (AtA[1][0]*AtA[2][2] - AtA[1][2]*AtA[2][0]) +
        AtA[0][2]*(AtA[1][0]*Atb[2]    - Atb[1]   *AtA[2][0]);

      // Slope at centroid of (i,j,k)
      ys = detAtA_y / detAtA;

    } else {

      amrex::Real du_yl = 2.0*(state(i,j  ,k,n) - state(i,j-1,k,n));
      amrex::Real du_yr = 2.0*(state(i,j+1,k,n) - state(i,j  ,k,n));
      amrex::Real du_yc = 0.5*(state(i,j+1,k,n) - state(i,j-1,k,n));

      amrex::Real yslope = amrex::min(std::abs(du_yl), std::abs(du_yc), std::abs(du_yr));
      yslope = (du_yr*du_yl > 0.0) ? yslope : 0.0;
      ys = (du_yc > 0.0) ? yslope : -yslope;
    }

    // Z direction
    if(flag(i,j,k  ).isSingleValued() or
      (flag(i,j,k-1).isSingleValued() or not flag(i,j,k).isConnected(0,0,-1)) or
      (flag(i,j,k+1).isSingleValued() or not flag(i,j,k).isConnected(0,0, 1))) 
    {

      amrex::Real detAtA_z =
        AtA[0][0]*(AtA[1][1]*Atb[2]    - Atb[1]   *AtA[1][2]) -
        AtA[0][1]*(AtA[1][0]*Atb[2]    - Atb[1]   *AtA[2][0]) +
        Atb[0]   *(AtA[1][0]*AtA[2][1] - AtA[1][1]*AtA[2][0]);

      zs = detAtA_z / detAtA;

    } else {

      amrex::Real du_zl = 2.0*(state(i,j,k  ,n) - state(i,j,k-1,n));
      amrex::Real du_zr = 2.0*(state(i,j,k+1,n) - state(i,j,k  ,n));
      amrex::Real du_zc = 0.5*(state(i,j,k+1,n) - state(i,j,k-1,n));

      amrex::Real zslope = amrex::min(std::abs(du_zl), std::abs(du_zc), std::abs(du_zr));
      zslope = (du_zr*du_zl > 0.0) ? zslope : 0.0;
      zs = (du_zc > 0.0) ? zslope : -zslope;
    }

   return {xs,ys,zs};
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::GpuArray<amrex::Real,3>
incflo_slopes_extdir_eb (int i, int j, int k, int n,
                  amrex::Array4<amrex::Real const> const& state,
                  amrex::Array4<amrex::Real const> const& ccent,
                  amrex::Array4<amrex::EBCellFlag const> const& flag,
                  bool edlo_x, bool edhi_x, int domlo_x, int domhi_x,
                  bool edlo_y, bool edhi_y, int domlo_y, int domhi_y,
                  bool edlo_z, bool edhi_z, int domlo_z, int domhi_z) noexcept
{
    // First get EB-aware slope that doesn't know about extdir
    amrex::GpuArray<amrex::Real,3> eb_slope = incflo_slopes_eb (i,j,k,n,state,ccent,flag);
    amrex::Real xslope = eb_slope[0];
    amrex::Real yslope = eb_slope[1];
    amrex::Real zslope = eb_slope[2];

    // Now correct only those cells which are not affected by EB but are by extdir
    if        (edlo_x and i == domlo_x) {
        if( !( flag(i  ,j,k).isSingleValued() or
              (flag(i-1,j,k).isSingleValued() or not flag(i,j,k).isConnected(-1,0,0)) or
              (flag(i+1,j,k).isSingleValued() or not flag(i,j,k).isConnected( 1,0,0))) ) 
        {
           amrex::Real dl = 2.0*(state(i  ,j,k,n) - state(i-1,j,k,n));
           amrex::Real dr = 2.0*(state(i+1,j,k,n) - state(i  ,j,k,n));
           amrex::Real dc = 
              (state(i+1,j,k,n)+3.0*state(i,j,k,n)-4.0*state(i-1,j,k,n))/3.0;

           amrex::Real slope = amrex::min(std::abs(dl),std::abs(dc),std::abs(dr));
           slope = (dr*dl > 0.0) ? slope : 0.0;
           xslope = (dc > 0.0) ? slope : -slope;
        }

    } else if (edhi_x and i == domhi_x) {
        if( !( flag(i  ,j,k).isSingleValued() or
              (flag(i-1,j,k).isSingleValued() or not flag(i,j,k).isConnected(-1,0,0)) or
              (flag(i+1,j,k).isSingleValued() or not flag(i,j,k).isConnected( 1,0,0))) )
        {
           amrex::Real dl = 2.0*(state(i  ,j,k,n) - state(i-1,j,k,n));
           amrex::Real dr = 2.0*(state(i+1,j,k,n) - state(i  ,j,k,n));
           amrex::Real dc = 
              (4.0*state(i+1,j,k,n)-3.0*state(i,j,k,n)-state(i-1,j,k,n))/3.0;

           amrex::Real slope = amrex::min(std::abs(dl),std::abs(dc),std::abs(dr));
           slope = (dr*dl > 0.0) ? slope : 0.0;
           xslope = (dc > 0.0) ? slope : -slope;
        }
    }

    if        (edlo_y and j == domlo_y) {
        if( !( flag(i,j  ,k).isSingleValued() or
              (flag(i,j-1,k).isSingleValued() or not flag(i,j,k).isConnected(0,-1,0)) or
              (flag(i,j+1,k).isSingleValued() or not flag(i,j,k).isConnected(0, 1,0))) )
        {
           amrex::Real dl = 2.0*(state(i,j  ,k,n) - state(i,j-1,k,n));
           amrex::Real dr = 2.0*(state(i,j+1,k,n) - state(i,j  ,k,n));
           amrex::Real dc = 
              (state(i,j+1,k,n)+3.0*state(i,j,k,n)-4.0*state(i,j-1,k,n))/3.0;

           amrex::Real slope = amrex::min(std::abs(dl),std::abs(dc),std::abs(dr));
           slope = (dr*dl > 0.0) ? slope : 0.0;
           yslope = (dc > 0.0) ? slope : -slope;
        }
    } else if (edhi_y and j == domhi_y) {
        if( !(flag(i,j  ,k).isSingleValued() or
           (flag(i,j-1,k).isSingleValued() or not flag(i,j,k).isConnected(0,-1,0)) or
           (flag(i,j+1,k).isSingleValued() or not flag(i,j,k).isConnected(0, 1,0))) )
        {
           amrex::Real dl = 2.0*(state(i,j  ,k,n) - state(i,j-1,k,n));
           amrex::Real dr = 2.0*(state(i,j+1,k,n) - state(i,j  ,k,n));
           amrex::Real dc = 
              (4.0*state(i,j+1,k,n)-3.0*state(i,j,k,n)-state(i,j-1,k,n))/3.0;

           amrex::Real slope = amrex::min(std::abs(dl),std::abs(dc),std::abs(dr));
           slope = (dr*dl > 0.0) ? slope : 0.0;
           yslope = (dc > 0.0) ? slope : -slope;
        }
    }

    if        (edlo_z and k == domlo_z) {
        if( !( flag(i,j,k  ).isSingleValued() or
              (flag(i,j,k-1).isSingleValued() or not flag(i,j,k).isConnected(0,0,-1)) or
              (flag(i,j,k+1).isSingleValued() or not flag(i,j,k).isConnected(0,0, 1))) )
        {
           amrex::Real dl = 2.0*(state(i,j,k  ,n) - state(i,j,k-1,n));
           amrex::Real dr = 2.0*(state(i,j,k+1,n) - state(i,j,k  ,n));
           amrex::Real dc = 
              (state(i,j,k+1,n)+3.0*state(i,j,k,n)-4.0*state(i,j,k-1,n))/3.0;

           amrex::Real slope = amrex::min(std::abs(dl),std::abs(dc),std::abs(dr));
           slope = (dr*dl > 0.0) ? slope : 0.0;
           zslope = (dc > 0.0) ? slope : -slope;
        }
    } else if (edhi_z and k == domhi_z) {
        if( !( flag(i,j,k  ).isSingleValued() or
             (flag(i,j,k-1).isSingleValued() or not flag(i,j,k).isConnected(0,0,-1)) or
             (flag(i,j,k+1).isSingleValued() or not flag(i,j,k).isConnected(0,0, 1))) )
        {
           amrex::Real dl = 2.0*(state(i,j,k  ,n) - state(i,j,k-1,n));
           amrex::Real dr = 2.0*(state(i,j,k+1,n) - state(i,j,k  ,n));
           amrex::Real dc = 
              (4.0*state(i,j,k+1,n)-3.0*state(i,j,k,n)-state(i,j,k-1,n))/3.0;

           amrex::Real slope = amrex::min(std::abs(dl),std::abs(dc),std::abs(dr));
           slope = (dr*dl > 0.0) ? slope : 0.0;
           zslope = (dc > 0.0) ? slope : -slope;
        }
    }

    return {xslope,yslope,zslope};
}

#endif

}
#endif
