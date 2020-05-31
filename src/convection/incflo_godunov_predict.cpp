#include "amr-wind/convection/incflo_godunov_plm.H"
#include "amr-wind/convection/incflo_godunov_ppm.H"

#include <AMReX_BCRec.H>
#include "amr-wind/convection/Godunov.H"

using namespace amrex;

void godunov::make_trans_velocities (int lev, Box const& xbx, Box const& ybx, Box const& zbx,
                                    Array4<Real> const& u_ad,
                                    Array4<Real> const& v_ad,
                                    Array4<Real> const& w_ad,
                                    Array4<Real const> const& Imx,
                                    Array4<Real const> const& Ipx,
                                    Array4<Real const> const& Imy,
                                    Array4<Real const> const& Ipy,
                                    Array4<Real const> const& Imz,
                                    Array4<Real const> const& Ipz,
                                    Array4<Real const> const& vel,
                                    Array4<Real const> const& f,
                                    Vector<Geometry> geom,
                                    Real dt,
                                    amrex::Gpu::DeviceVector<amrex::BCRec>& bcrec_device,
                                    bool godunov_use_forces_in_trans)

{
    BL_PROFILE("amr-wind::godunov::make_trans_velocities")
    Real l_dt = dt;
    bool l_use_forces_in_trans = godunov_use_forces_in_trans;

    const Box& domain = geom[lev].Domain();
    const Dim3 dlo = amrex::lbound(domain);
    const Dim3 dhi = amrex::ubound(domain);

    BCRec const* pbc = bcrec_device.data();

    amrex::ParallelFor(xbx, ybx, zbx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        // We only care about x-velocity on x-faces here
        constexpr int n = 0;

        Real lo, hi;
        if (l_use_forces_in_trans) {
            lo = Ipx(i-1,j,k,n) + 0.5*l_dt*f(i-1,j,k,n);
            hi = Imx(i  ,j,k,n) + 0.5*l_dt*f(i  ,j,k,n);
        } else {
            lo = Ipx(i-1,j,k,n);
            hi = Imx(i  ,j,k,n);
        }

        auto bc = pbc[n];
        Godunov_trans_xbc(i, j, k, n, vel, lo, hi, lo, bc.lo(0), bc.hi(0), dlo.x, dhi.x);

        constexpr Real small_vel = 1e-10;

        Real st = ( (lo+hi) >= 0.) ? lo : hi;
        bool ltm = ( (lo <= 0. && hi >= 0.) || (amrex::Math::abs(lo+hi) < small_vel) );
        u_ad(i,j,k) = ltm ? 0. : st;
    },
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        // We only care about y-velocity on y-faces here
        constexpr int n = 1;

        Real lo, hi;
        if (l_use_forces_in_trans) {
            lo = Ipy(i,j-1,k,n) + 0.5*l_dt*f(i,j-1,k,n);
            hi = Imy(i,j  ,k,n) + 0.5*l_dt*f(i,j  ,k,n);
        } else {
            lo = Ipy(i,j-1,k,n);
            hi = Imy(i,j  ,k,n);
        }

        auto bc = pbc[n];
        Godunov_trans_ybc(i, j, k, n, vel, lo, hi, lo, bc.lo(1), bc.hi(1), dlo.y, dhi.y);

        constexpr Real small_vel = 1e-10;

        Real st = ( (lo+hi) >= 0.) ? lo : hi;
        bool ltm = ( (lo <= 0. && hi >= 0.) || (amrex::Math::abs(lo+hi) < small_vel) );
        v_ad(i,j,k) = ltm ? 0. : st;
    },
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        // We only care about z-velocity on z-faces here
        constexpr int n = 2;

        Real lo, hi;
        if (l_use_forces_in_trans) {
            lo = Ipz(i,j,k-1,n) + 0.5*l_dt*f(i,j,k-1,n);
            hi = Imz(i,j,k  ,n) + 0.5*l_dt*f(i,j,k  ,n);
        } else {
            lo = Ipz(i,j,k-1,n);
            hi = Imz(i,j,k  ,n);
        }

        auto bc = pbc[n];
        Godunov_trans_zbc(i, j, k, n, vel, lo, hi, lo, bc.lo(2), bc.hi(2), dlo.z, dhi.z);

        constexpr Real small_vel = 1e-10;

        Real st = ( (lo+hi) >= 0.) ? lo : hi;
        bool ltm = ( (lo <= 0. && hi >= 0.) || (amrex::Math::abs(lo+hi) < small_vel) );
        w_ad(i,j,k) = ltm ? 0. : st;
    });
}

void godunov::predict_godunov (int lev, Box const& bx, int ncomp,
                              Box const& xbx, Box const& ybx, Box const& zbx,
                              Array4<Real> const& qx,
                              Array4<Real> const& qy,
                              Array4<Real> const& qz,
                              Array4<Real const> const& q,
                              Array4<Real const> const& u_ad,
                              Array4<Real const> const& v_ad,
                              Array4<Real const> const& w_ad,
                              Array4<Real> const& Imx,
                              Array4<Real> const& Ipx,
                              Array4<Real> const& Imy,
                              Array4<Real> const& Ipy,
                              Array4<Real> const& Imz,
                              Array4<Real> const& Ipz,
                              Array4<Real const> const& f,
                              Real* p,
                              Vector<Geometry> geom,
                              Real dt,
                              amrex::Gpu::DeviceVector<amrex::BCRec>& bcrec_device,
                              bool godunov_use_forces_in_trans)
{
    BL_PROFILE("amr-wind::godunov::predict_godunov")
    Real l_dt = dt;
    bool l_use_forces_in_trans = godunov_use_forces_in_trans;

    const Box& domain = geom[lev].Domain();
    const Dim3 dlo = amrex::lbound(domain);
    const Dim3 dhi = amrex::ubound(domain);
    Real dx = geom[lev].CellSize(0);
    Real dy = geom[lev].CellSize(1);
    Real dz = geom[lev].CellSize(2);

    BCRec const* pbc = bcrec_device.data();

    Box xebox = Box(bx).grow(1,1).grow(2,1).surroundingNodes(0);
    Box yebox = Box(bx).grow(0,1).grow(2,1).surroundingNodes(1);
    Box zebox = Box(bx).grow(0,1).grow(1,1).surroundingNodes(2);
    Array4<Real> xlo = makeArray4(p, xebox, ncomp);
    p += xlo.size();
    Array4<Real> xhi = makeArray4(p, xebox, ncomp);
    p += xhi.size();
    Array4<Real> ylo = makeArray4(p, yebox, ncomp);
    p += ylo.size();
    Array4<Real> yhi = makeArray4(p, yebox, ncomp);
    p += yhi.size();
    Array4<Real> zlo = makeArray4(p, zebox, ncomp);
    p += zlo.size();
    Array4<Real> zhi = makeArray4(p, zebox, ncomp);
    p += zhi.size();

    amrex::ParallelFor(
        xebox, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            Real lo, hi;
            if (l_use_forces_in_trans) {
                lo = Ipx(i-1,j,k,n) + 0.5*l_dt*f(i-1,j,k,n);
                hi = Imx(i  ,j,k,n) + 0.5*l_dt*f(i  ,j,k,n);
            } else {
                lo = Ipx(i-1,j,k,n);
                hi = Imx(i  ,j,k,n);
            }

            Real uad = u_ad(i,j,k);
            auto bc = pbc[n];

            Godunov_trans_xbc(i, j, k, n, q, lo, hi, uad, bc.lo(0), bc.hi(0), dlo.x, dhi.x);

            xlo(i,j,k,n) = lo;
            xhi(i,j,k,n) = hi;

            constexpr Real small_vel = 1e-10;

            Real st = (uad >= 0.) ? lo : hi;
            Real fu = (amrex::Math::abs(uad) < small_vel) ? 0.0 : 1.0;
            Imx(i, j, k, n) = fu*st + (1.0 - fu) *0.5 * (hi + lo); // store xedge
        },
        yebox, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            Real lo, hi;
            if (l_use_forces_in_trans) {
                lo = Ipy(i,j-1,k,n) + 0.5*l_dt*f(i,j-1,k,n);
                hi = Imy(i,j  ,k,n) + 0.5*l_dt*f(i,j  ,k,n);
            } else {
                lo = Ipy(i,j-1,k,n);
                hi = Imy(i,j  ,k,n);
            }

            Real vad = v_ad(i,j,k);
            auto bc = pbc[n];

            Godunov_trans_ybc(i, j, k, n, q, lo, hi, vad, bc.lo(1), bc.hi(1), dlo.y, dhi.y);

            ylo(i,j,k,n) = lo;
            yhi(i,j,k,n) = hi;

            constexpr Real small_vel = 1e-10;

            Real st = (vad >= 0.) ? lo : hi;
            Real fu = (amrex::Math::abs(vad) < small_vel) ? 0.0 : 1.0;
            Imy(i, j, k, n) = fu*st + (1.0 - fu)*0.5*(hi + lo); // store yedge
        },
        zebox, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            Real lo, hi;
            if (l_use_forces_in_trans) {
                lo = Ipz(i,j,k-1,n) + 0.5*l_dt*f(i,j,k-1,n);
                hi = Imz(i,j,k  ,n) + 0.5*l_dt*f(i,j,k  ,n);
            } else {
                lo = Ipz(i,j,k-1,n);
                hi = Imz(i,j,k  ,n);
            }

            Real wad = w_ad(i,j,k);
            auto bc = pbc[n];
            
            Godunov_trans_zbc(i, j, k, n, q, lo, hi, wad, bc.lo(2), bc.hi(2), dlo.z, dhi.z);

            zlo(i,j,k,n) = lo;
            zhi(i,j,k,n) = hi;

            constexpr Real small_vel = 1e-10;

            Real st = (wad >= 0.) ? lo : hi;
            Real fu = (amrex::Math::abs(wad) < small_vel) ? 0.0 : 1.0;
            Imz(i, j, k, n) = fu*st + (1.0 - fu)*0.5*(hi + lo); // store zedge
        });

    Array4<Real> xedge = Imx;
    Array4<Real> yedge = Imy;
    Array4<Real> zedge = Imz;

    Array4<Real> divu = makeArray4(Ipx.dataPtr(), grow(bx,1), 1);
    amrex::ParallelFor(Box(divu), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        divu(i,j,k) = 0.0;
    });

    // We can reuse the space in Ipy and Ipz.

    //
    // X-Flux
    //
    Box const xbxtmp = Box(xbx).enclosedCells().grow(0,1);
    Array4<Real> yzlo = makeArray4(Ipy.dataPtr(), amrex::surroundingNodes(xbxtmp,1), 1);
    Array4<Real> zylo = makeArray4(Ipz.dataPtr(), amrex::surroundingNodes(xbxtmp,2), 1);
    // Add d/dy term to z-faces
    // Start with {zlo,zhi} --> {zylo, zyhi} and upwind using w_ad to {zylo}
    // Add d/dz to y-faces
    // Start with {ylo,yhi} --> {yzlo, yzhi} and upwind using v_ad to {yzlo}
    amrex::ParallelFor(Box(zylo), Box(yzlo),
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        constexpr int n = 0;
        const auto bc = pbc[n]; 
        Real l_zylo, l_zyhi;
        Godunov_corner_couple_zy(l_zylo, l_zyhi,
                                 i, j, k, n, l_dt, dy, false,
                                 zlo(i,j,k,n), zhi(i,j,k,n),
                                 q, divu, v_ad, yedge);

        Real wad = w_ad(i,j,k);
        Godunov_trans_zbc(i, j, k, n, q, l_zylo, l_zyhi, wad, bc.lo(2), bc.hi(2), dlo.z, dhi.z);

        constexpr Real small_vel = 1.e-10;

        Real st = (wad >= 0.) ? l_zylo : l_zyhi;
        Real fu = (amrex::Math::abs(wad) < small_vel) ? 0.0 : 1.0;
        zylo(i,j,k) = fu*st + (1.0 - fu) * 0.5 * (l_zyhi + l_zylo);
    },
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        constexpr int n = 0;
        const auto bc = pbc[n];
        Real l_yzlo, l_yzhi;
        Godunov_corner_couple_yz(l_yzlo, l_yzhi,
                                 i, j, k, n, l_dt, dz, false,
                                 ylo(i,j,k,n), yhi(i,j,k,n),
                                 q, divu, w_ad, zedge);

        Real vad = v_ad(i,j,k);
        Godunov_trans_ybc(i, j, k, n, q, l_yzlo, l_yzhi, vad, bc.lo(1), bc.hi(1), dlo.y, dhi.y);

        constexpr Real small_vel = 1.e-10;

        Real st = (vad >= 0.) ? l_yzlo : l_yzhi;
        Real fu = (amrex::Math::abs(vad) < small_vel) ? 0.0 : 1.0;
        yzlo(i,j,k) = fu*st + (1.0 - fu) * 0.5 * (l_yzhi + l_yzlo);
    });
    //
    amrex::ParallelFor(xbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        constexpr int n = 0;
        auto bc = pbc[n];
        Real stl = xlo(i,j,k,n) - (0.25*l_dt/dy)*(v_ad(i-1,j+1,k  )+v_ad(i-1,j,k))*
                                                 (yzlo(i-1,j+1,k  )-yzlo(i-1,j,k))
                                - (0.25*l_dt/dz)*(w_ad(i-1,j  ,k+1)+w_ad(i-1,j,k))*
                                                 (zylo(i-1,j  ,k+1)-zylo(i-1,j,k));
        Real sth = xhi(i,j,k,n) - (0.25*l_dt/dy)*(v_ad(i  ,j+1,k  )+v_ad(i  ,j,k))*
                                                 (yzlo(i  ,j+1,k  )-yzlo(i  ,j,k))
                                - (0.25*l_dt/dz)*(w_ad(i  ,j  ,k+1)+w_ad(i  ,j,k))*
                                                 (zylo(i  ,j  ,k+1)-zylo(i  ,j,k));
        if (!l_use_forces_in_trans) {
            stl += 0.5 * l_dt * f(i-1,j,k,n);
            sth += 0.5 * l_dt * f(i  ,j,k,n);
        }

        Godunov_cc_xbc_lo(i, j, k, n, q, stl, sth, u_ad, bc.lo(0), dlo.x);
        Godunov_cc_xbc_hi(i, j, k, n, q, stl, sth, u_ad, bc.hi(0), dhi.x);

        constexpr Real small_vel = 1.e-10;

        Real st = ( (stl+sth) >= 0.) ? stl : sth;
        bool ltm = ( (stl <= 0. && sth >= 0.) || (amrex::Math::abs(stl+sth) < small_vel) );
        qx(i,j,k) = ltm ? 0. : st;
    });

    //
    // Y-Flux
    //
    Box const ybxtmp = Box(ybx).enclosedCells().grow(1,1);
    Array4<Real> xzlo = makeArray4(Ipy.dataPtr(), amrex::surroundingNodes(ybxtmp,0), 1);
    Array4<Real> zxlo = makeArray4(Ipz.dataPtr(), amrex::surroundingNodes(ybxtmp,2), 1);

    // Add d/dz to x-faces
    // Start with {xlo,xhi} --> {xzlo, xzhi} and upwind using u_ad to {xzlo}
    // Add d/dx term to z-faces
    // Start with {zlo,zhi} --> {zxlo, zxhi} and upwind using w_ad to {zxlo}
    amrex::ParallelFor(Box(xzlo), Box(zxlo),
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        constexpr int n = 1;
        const auto bc = pbc[n];
        Real l_xzlo, l_xzhi;
        Godunov_corner_couple_xz(l_xzlo, l_xzhi,
                                 i, j, k, n, l_dt, dz, false,
                                 xlo(i,j,k,n),  xhi(i,j,k,n),
                                 q, divu, w_ad, zedge);

        Real uad = u_ad(i,j,k);
        Godunov_trans_xbc(i, j, k, n, q, l_xzlo, l_xzhi, uad, bc.lo(0), bc.hi(0), dlo.x, dhi.x);

        constexpr Real small_vel = 1.e-10;

        Real st = (uad >= 0.) ? l_xzlo : l_xzhi;
        Real fu = (amrex::Math::abs(uad) < small_vel) ? 0.0 : 1.0;
        xzlo(i,j,k) = fu*st + (1.0 - fu) * 0.5 * (l_xzhi + l_xzlo);
    },
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        constexpr int n = 1;
        const auto bc = pbc[n];
        Real l_zxlo, l_zxhi;
        Godunov_corner_couple_zx(l_zxlo, l_zxhi,
                                 i, j, k, n, l_dt, dx, false,
                                 zlo(i,j,k,n), zhi(i,j,k,n),
                                 q, divu, u_ad, xedge);

        Real wad = w_ad(i,j,k);
        Godunov_trans_zbc(i, j, k, n, q, l_zxlo, l_zxhi, wad, bc.lo(2), bc.hi(2), dlo.z, dhi.z);

        constexpr Real small_vel = 1.e-10;

        Real st = (wad >= 0.) ? l_zxlo : l_zxhi;
        Real fu = (amrex::Math::abs(wad) < small_vel) ? 0.0 : 1.0;
        zxlo(i,j,k) = fu*st + (1.0 - fu) * 0.5 * (l_zxhi + l_zxlo);
    });
    //
    amrex::ParallelFor(ybx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        constexpr int n = 1;
        auto bc = pbc[n];
        Real stl = ylo(i,j,k,n) - (0.25*l_dt/dx)*(u_ad(i+1,j-1,k  )+u_ad(i,j-1,k))*
                                                 (xzlo(i+1,j-1,k  )-xzlo(i,j-1,k))
                                - (0.25*l_dt/dz)*(w_ad(i  ,j-1,k+1)+w_ad(i,j-1,k))*
                                                 (zxlo(i  ,j-1,k+1)-zxlo(i,j-1,k));
        Real sth = yhi(i,j,k,n) - (0.25*l_dt/dx)*(u_ad(i+1,j  ,k  )+u_ad(i,j  ,k))*
                                                 (xzlo(i+1,j  ,k  )-xzlo(i,j  ,k))
                                - (0.25*l_dt/dz)*(w_ad(i  ,j  ,k+1)+w_ad(i,j  ,k))*
                                                 (zxlo(i  ,j  ,k+1)-zxlo(i,j  ,k));
        if (!l_use_forces_in_trans) {
           stl += 0.5 * l_dt * f(i,j-1,k,n);
           sth += 0.5 * l_dt * f(i,j  ,k,n);
        }

        Godunov_cc_ybc_lo(i, j, k, n, q, stl, sth, v_ad, bc.lo(1), dlo.y);
        Godunov_cc_ybc_hi(i, j, k, n, q, stl, sth, v_ad, bc.hi(1), dhi.y);

        constexpr Real small_vel = 1.e-10;

        Real st = ( (stl+sth) >= 0.) ? stl : sth;
        bool ltm = ( (stl <= 0. && sth >= 0.) || (amrex::Math::abs(stl+sth) < small_vel) );
        qy(i,j,k) = ltm ? 0. : st;
    });

    //
    // Z-Flux
    //
    Box const zbxtmp = Box(zbx).enclosedCells().grow(2,1);
    Array4<Real> xylo = makeArray4(Ipy.dataPtr(), amrex::surroundingNodes(zbxtmp,0), 1);
    Array4<Real> yxlo = makeArray4(Ipz.dataPtr(), amrex::surroundingNodes(zbxtmp,1), 1);

    amrex::ParallelFor(Box(xylo), Box(yxlo),
    //
    // Add d/dy term to x-faces
    // Start with {xlo,xhi} --> {xylo, xyhi} and upwind using u_ad to {xylo}
    //
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        constexpr int n = 2;
        const auto bc = pbc[n];
        Real l_xylo, l_xyhi;
        Godunov_corner_couple_xy(l_xylo, l_xyhi,
                                 i, j, k, n, l_dt, dy, false,
                                 xlo(i,j,k,n), xhi(i,j,k,n),
                                 q, divu, v_ad, yedge);

        Real uad = u_ad(i,j,k);
        Godunov_trans_xbc(i, j, k, n, q, l_xylo, l_xyhi, uad, bc.lo(0), bc.hi(0), dlo.x, dhi.x);

        constexpr Real small_vel = 1.e-10;

        Real st = (uad >= 0.) ? l_xylo : l_xyhi;
        Real fu = (amrex::Math::abs(uad) < small_vel) ? 0.0 : 1.0;
        xylo(i,j,k) = fu*st + (1.0 - fu) * 0.5 * (l_xyhi + l_xylo);
    },
    //
    // Add d/dx term to y-faces
    // Start with {ylo,yhi} --> {yxlo, yxhi} and upwind using v_ad to {yxlo}
    //
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        constexpr int n = 2;
        const auto bc = pbc[n];
        Real l_yxlo, l_yxhi;
        Godunov_corner_couple_yx(l_yxlo, l_yxhi,
                                 i, j, k, n, l_dt, dx, false,
                                 ylo(i,j,k,n), yhi(i,j,k,n),
                                 q, divu, u_ad, xedge);

        Real vad = v_ad(i,j,k);
        Godunov_trans_ybc(i, j, k, n, q, l_yxlo, l_yxhi, vad, bc.lo(1), bc.hi(1), dlo.y, dhi.y);

        constexpr Real small_vel = 1.e-10;

        Real st = (vad >= 0.) ? l_yxlo : l_yxhi;
        Real fu = (amrex::Math::abs(vad) < small_vel) ? 0.0 : 1.0;
        yxlo(i,j,k) = fu*st + (1.0 - fu) * 0.5 * (l_yxhi + l_yxlo);
    });
    //
    amrex::ParallelFor(zbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        constexpr int n = 2;
        auto bc = pbc[n];
        Real stl = zlo(i,j,k,n) - (0.25*l_dt/dx)*(u_ad(i+1,j  ,k-1)+u_ad(i,j,k-1))*
                                                 (xylo(i+1,j  ,k-1)-xylo(i,j,k-1))
                                - (0.25*l_dt/dy)*(v_ad(i  ,j+1,k-1)+v_ad(i,j,k-1))*
                                                 (yxlo(i  ,j+1,k-1)-yxlo(i,j,k-1));
        Real sth = zhi(i,j,k,n) - (0.25*l_dt/dx)*(u_ad(i+1,j  ,k  )+u_ad(i,j,k  ))*
                                                 (xylo(i+1,j  ,k  )-xylo(i,j,k  ))
                                - (0.25*l_dt/dy)*(v_ad(i  ,j+1,k  )+v_ad(i,j,k  ))*
                                                 (yxlo(i  ,j+1,k  )-yxlo(i,j,k  ));

        if (!l_use_forces_in_trans) {
           stl += 0.5 * l_dt * f(i,j,k-1,n);
           sth += 0.5 * l_dt * f(i,j,k  ,n);
        }

        Godunov_cc_zbc_lo(i, j, k, n, q, stl, sth, w_ad, bc.lo(2), dlo.z);
        Godunov_cc_zbc_hi(i, j, k, n, q, stl, sth, w_ad, bc.hi(2), dhi.z);

        constexpr Real small_vel = 1.e-10;

        Real st = ( (stl+sth) >= 0.) ? stl : sth;
        bool ltm = ( (stl <= 0. && sth >= 0.) || (amrex::Math::abs(stl+sth) < small_vel) );
        qz(i,j,k) = ltm ? 0. : st;
    });
}
