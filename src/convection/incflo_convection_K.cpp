#include <incflo_convection_K.H>

using namespace amrex;

void incflo_predict_vels_on_faces (Box const& ubx, Box const& vbx, Box const& wbx,
                                   Array4<Real> const& u, Array4<Real> const& v,
                                   Array4<Real> const& w, Array4<Real const> const& vcc)
{
    constexpr Real small = 1.e-10;

    amrex::ParallelFor(ubx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        Real upls = vcc(i  ,j,k,0) - 0.5 * incflo_xslope(i  ,j,k,0,vcc);
        Real umns = vcc(i-1,j,k,0) + 0.5 * incflo_xslope(i-1,j,k,0,vcc);
        if (umns < 0.0 and upls > 0.0) {
            u(i,j,k) = 0.0;
        } else {
            Real avg = 0.5 * (upls + umns);
            if (std::abs(avg) < small) {
                u(i,j,k) = 0.0;
            } else if (avg > 0.0) {
                u(i,j,k) = umns;
            } else {
                u(i,j,k) = upls;
            }
        }
    });

    amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        Real vpls = vcc(i  ,j,k,1) - 0.5 * incflo_xslope(i  ,j,k,1,vcc);
        Real vmns = vcc(i-1,j,k,1) + 0.5 * incflo_xslope(i-1,j,k,1,vcc);
        if (vmns < 0.0 and vpls > 0.0) {
            v(i,j,k) = 0.0;
        } else {
            Real avg = 0.5 * (vpls + vmns);
            if (std::abs(avg) < small) {
                v(i,j,k) = 0.0;
            } else if (avg > 0.0) {
                v(i,j,k) = vmns;
            } else {
                v(i,j,k) = vpls;
            }
        }
    });

    amrex::ParallelFor(wbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        Real wpls = vcc(i  ,j,k,2) - 0.5 * incflo_xslope(i  ,j,k,2,vcc);
        Real wmns = vcc(i-1,j,k,2) + 0.5 * incflo_xslope(i-1,j,k,2,vcc);
        if (wmns < 0.0 and wpls > 0.0) {
            w(i,j,k) = 0.0;
        } else {
            Real avg = 0.5 * (wpls + wmns);
            if (std::abs(avg) < small) {
                w(i,j,k) = 0.0;
            } else if (avg > 0.0) {
                w(i,j,k) = wmns;
            } else {
                w(i,j,k) = wpls;
            }
        }
    });
}
