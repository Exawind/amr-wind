#include <incflo_godunov_ppm.H>
#include <incflo.H>

using namespace amrex;

void incflo::predict_ppm (int lev, Box const& bx, int ncomp,
                          Array4<Real> const& Imx,
                          Array4<Real> const& Ipx,
                          Array4<Real> const& Imy,
                          Array4<Real> const& Ipy,
                          Array4<Real> const& Imz,
                          Array4<Real> const& Ipz,
                          Array4<Real const> const& q,
                          Array4<Real const> const& vel)
{
    Real l_dt = m_dt;
    const auto dx = Geom(lev).CellSizeArray();
    const Box& domain = Geom(lev).Domain();
    const Dim3 dlo = amrex::lbound(domain);
    const Dim3 dhi = amrex::ubound(domain);

    BCRec const* pbc = get_velocity_bcrec_device_ptr();

    amrex::ParallelFor(bx, AMREX_SPACEDIM,
    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        Godunov_ppm_pred_x(i,j,k,n,l_dt,dx[0],q,vel,Imx,Ipx,pbc[0],dlo.x,dhi.x);
        Godunov_ppm_pred_y(i,j,k,n,l_dt,dx[1],q,vel,Imy,Ipy,pbc[1],dlo.y,dhi.y);
        Godunov_ppm_pred_z(i,j,k,n,l_dt,dx[2],q,vel,Imz,Ipz,pbc[2],dlo.z,dhi.z);
    });
}
