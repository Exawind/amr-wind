#include <incflo_godunov_ppm.H>
#include <incflo.H>

using namespace amrex;

void incflo::predict_ppm (int lev, Box const& bx, int /* ncomp */,
                          Array4<Real> const& Imx,
                          Array4<Real> const& Ipx,
                          Array4<Real> const& Imy,
                          Array4<Real> const& Ipy,
                          Array4<Real> const& Imz,
                          Array4<Real> const& Ipz,
                          Array4<Real const> const& q,
                          Array4<Real const> const& vel)
{
    const auto dx = Geom(lev).CellSizeArray();
    const Box& domain = Geom(lev).Domain();
    const Dim3 dlo = amrex::lbound(domain);
    const Dim3 dhi = amrex::ubound(domain);

    Real l_dtdx = m_time.deltaT() / dx[0];
    Real l_dtdy = m_time.deltaT() / dx[1];
    Real l_dtdz = m_time.deltaT() / dx[2];

    BCRec const* pbc = velocity().bcrec_device().data();

    amrex::ParallelFor(bx, AMREX_SPACEDIM, 
    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        Godunov_ppm_pred_x(i,j,k,n,l_dtdx,vel(i,j,k,0),q,Imx,Ipx,pbc[n],dlo.x,dhi.x);
        Godunov_ppm_pred_y(i,j,k,n,l_dtdy,vel(i,j,k,1),q,Imy,Ipy,pbc[n],dlo.y,dhi.y);
        Godunov_ppm_pred_z(i,j,k,n,l_dtdz,vel(i,j,k,2),q,Imz,Ipz,pbc[n],dlo.z,dhi.z);
    });
}
