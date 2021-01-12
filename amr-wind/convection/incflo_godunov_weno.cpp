#include "amr-wind/convection/incflo_godunov_weno.H"
#include "amr-wind/convection/Godunov.H"

using namespace amrex;

void godunov::predict_weno(
    int lev,
    Box const& bx,
    int /* ncomp */,
    Array4<Real> const& Imx,
    Array4<Real> const& Ipx,
    Array4<Real> const& Imy,
    Array4<Real> const& Ipy,
    Array4<Real> const& Imz,
    Array4<Real> const& Ipz,
    Array4<Real const> const& q,
    Array4<Real const> const& vel,
    Vector<Geometry> geom,
    Real dt,
    amrex::Gpu::DeviceVector<amrex::BCRec>& bcrec_device,
    bool weno_js)
{
    BL_PROFILE("amr-wind::godunov::predict_weno");
    const auto dx = geom[lev].CellSizeArray();
    const Box& domain = geom[lev].Domain();
    const Dim3 dlo = amrex::lbound(domain);
    const Dim3 dhi = amrex::ubound(domain);

    Real l_dtdx = dt / dx[0];
    Real l_dtdy = dt / dx[1];
    Real l_dtdz = dt / dx[2];

    BCRec const* pbc = bcrec_device.data();

    amrex::ParallelFor(
        bx, AMREX_SPACEDIM,
        [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
            Godunov_weno_pred_x(
                i, j, k, n, l_dtdx, vel(i, j, k, 0), q, Imx, Ipx, pbc[n], dlo.x,
                dhi.x, weno_js);
            Godunov_weno_pred_y(
                i, j, k, n, l_dtdy, vel(i, j, k, 1), q, Imy, Ipy, pbc[n], dlo.y,
                dhi.y, weno_js);
            Godunov_weno_pred_z(
                i, j, k, n, l_dtdz, vel(i, j, k, 2), q, Imz, Ipz, pbc[n], dlo.z,
                dhi.z, weno_js);
        });
}
