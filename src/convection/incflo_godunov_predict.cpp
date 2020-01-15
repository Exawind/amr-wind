#include "incflo_godunov_ppm.H" 
#include "incflo.H"
#include <AMReX_BCRec.H>
#include <iomanip> 

using namespace amrex;

void incflo::make_ppm_integrals (int lev, Box const& bx, int ncomp,
                                 Array4<Real> const& Imx,
                                 Array4<Real> const& Ipx,
                                 Array4<Real> const& Imy,
                                 Array4<Real> const& Ipy,
                                 Array4<Real> const& Imz,
                                 Array4<Real> const& Ipz,
                                 Array4<Real const> const& q,
                                 Array4<Real const> const& vel)
{
    Real l_dt = this->dt;
    bool l_use_forces_in_trans = this->use_forces_in_trans;

    const auto dx = Geom(lev).CellSize();
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

void incflo::make_trans_velocities (int lev, Box const& xbx, Box const& ybx, Box const& zbx,
                                    Array4<Real> const& u_ad,
                                    Array4<Real> const& v_ad,
                                    Array4<Real> const& w_ad,
                                    Array4<Real const> const& Imx,
                                    Array4<Real const> const& Ipx,
                                    Array4<Real const> const& Imy,
                                    Array4<Real const> const& Ipy,
                                    Array4<Real const> const& Imz,
                                    Array4<Real const> const& Ipz,
                                    Array4<Real const> const& vel)
{
    Real l_dt = this->dt;
    bool l_use_forces_in_trans = this->use_forces_in_trans;

    const Box& domain = Geom(lev).Domain();

    BCRec const* pbc = get_velocity_bcrec_device_ptr();

    amrex::ParallelFor(xbx, ybx, zbx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        // We only care about x-velocity on x-faces here
        constexpr int n = 0;

        Real lo, hi;
        if (l_use_forces_in_trans) {
            amrex::Abort("TODO: use_forces_in_trans");
            // lo = Ipx(i-1,j,k,n) + 0.5*l_dt*tf(i-1,j,k,n);
            // hi = Imx(i  ,j,k,n) + 0.5*l_dt*tf(i  ,j,k,n);
        } else {
            lo = Ipx(i-1,j,k,n);
            hi = Imx(i  ,j,k,n);
        }

        auto bc = pbc[n];
        Godunov_trans_xbc_lo(i, j, k, n, vel, lo, hi, lo, bc.lo(0),
                             domain.loVect(), domain.hiVect(), false, false);
        Godunov_trans_xbc_hi(i, j, k, n, vel, lo, hi, lo, bc.hi(0), 
                             domain.loVect(), domain.hiVect(), false, false);
        constexpr Real eps = 1e-6;

        Real st = ( (lo+hi) >= 0.) ? lo : hi;
        bool ltm = ( (lo <= 0. && hi >= 0.) || (std::abs(lo+hi) < eps) );
        u_ad(i,j,k) = ltm ? 0. : st;
    },
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        // We only care about y-velocity on y-faces here
        constexpr int n = 1;

        Real lo, hi;
        if (l_use_forces_in_trans) {
            amrex::Abort("TODO: use_forces_in_trans");
            // lo = Ipy(i,j-1,k,n) + 0.5*l_dt*tf(i,j-1,k,n);
            // hi = Imy(i,j,k,n)   + 0.5*l_dt*tf(i,j-1,k,n);
        } else {
            lo = Ipy(i,j-1,k,n);
            hi = Imy(i,j  ,k,n);
        }

        auto bc = pbc[n];
        Godunov_trans_ybc_lo(i, j, k, n, vel, lo, hi, lo, bc.lo(1),
                             domain.loVect(), domain.hiVect(), false, false);
        Godunov_trans_ybc_hi(i, j, k, n, vel, lo, hi, lo, bc.hi(1),
                             domain.loVect(), domain.hiVect(), false, false);
        constexpr Real eps = 1e-6;

        Real st = ( (lo+hi) >= 0.) ? lo : hi;
        bool ltm = ( (lo <= 0. && hi >= 0.) || (std::abs(lo+hi) < eps) );
        v_ad(i,j,k) = ltm ? 0. : st;
    },
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        // We only care about z-velocity on z-faces here
        constexpr int n = 2;

        Real lo, hi;
        if (l_use_forces_in_trans) {
            amrex::Abort("TODO: use_forces_in_trans");
            // lo = Ipz(i,j,k-1,n) + 0.5*l_dt*tf(i,j,k-1,n);
            // hi = Imz(i,j,k,n)   + 0.5*l_dt*tf(i,j,k-1,n);
        } else {
            lo = Ipz(i,j,k-1,n);
            hi = Imz(i,j,k  ,n);
        }

        auto bc = pbc[n];
        Godunov_trans_zbc_lo(i, j, k, n, vel, lo, hi, lo, bc.lo(2),
                             domain.loVect(), domain.hiVect(), false, false);
        Godunov_trans_zbc_hi(i, j, k, n, vel, lo, hi, lo, bc.hi(2),
                             domain.loVect(), domain.hiVect(), false, false);

        constexpr Real eps = 1e-6;

        Real st = ( (lo+hi) >= 0.) ? lo : hi;
        bool ltm = ( (lo <= 0. && hi >= 0.) || (std::abs(lo+hi) < eps) );
        w_ad(i,j,k) = ltm ? 0. : st;
    });
}

void incflo::predict_godunov_on_box (int lev, Box const& bx, int ncomp,
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
                                     Real* p)
{
    Real l_dt = this->dt;
    bool l_use_forces_in_trans = this->use_forces_in_trans;

    const Box& domain = Geom(lev).Domain();
    const Dim3 dlo = amrex::lbound(domain);
    const Dim3 dhi = amrex::ubound(domain);
    Real dx = Geom(lev).CellSize(0);
    Real dy = Geom(lev).CellSize(1);
    Real dz = Geom(lev).CellSize(2);

    // We only predict velocity in this routine, which is convectively not conservatively, updated
    GpuArray<int,3> iconserv{0,0,0};

    BCRec const* pbc = get_velocity_bcrec_device_ptr();

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
                // xxxxx lo = Ipx(i-1,j,k,n) + 0.5*l_dt*tf(i-1,j,k,n);
                // hi = Imx(i  ,j,k,n) + 0.5*l_dt*tf(i  ,j,k,n);
            } else {
                lo = Ipx(i-1,j,k,n);
                hi = Imx(i  ,j,k,n);
            }

            Real uad = u_ad(i,j,k);
            auto bc = pbc[n];

            Godunov_trans_xbc_lo(i, j, k, n, q, lo, hi, uad, bc.lo(0),
                                 domain.loVect(), domain.hiVect(), false, false);
            Godunov_trans_xbc_hi(i, j, k, n, q, lo, hi, uad, bc.hi(0),
                                 domain.loVect(), domain.hiVect(), false, false);
            xlo(i,j,k,n) = lo;
            xhi(i,j,k,n) = hi;

            constexpr Real eps = 1e-6;

            Real st = (uad >= 0.) ? lo : hi;
            Real fu = (std::abs(uad) < eps) ? 0.0 : 1.0;
            Imx(i, j, k, n) = fu*st + (1.0 - fu) *0.5 * (hi + lo); // store xedge
        },
        yebox, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            Real lo, hi;
            if (l_use_forces_in_trans) {
                // xxxxx lo = Ipy(i,j-1,k,n) + 0.5*l_dt*tf(i,j-1,k,n);
                // hi = Imy(i,j,k,n)   + 0.5*l_dt*tf(i,j-1,k,n);
            } else {
                lo = Ipy(i,j-1,k,n);
                hi = Imy(i,j  ,k,n);
            }

            Real vad = v_ad(i,j,k);
            auto bc = pbc[n];

            Godunov_trans_ybc_lo(i, j, k, n, q, lo, hi, vad, bc.lo(1),
                                 domain.loVect(), domain.hiVect(), false, false);
            Godunov_trans_ybc_hi(i, j, k, n, q, lo, hi, vad, bc.hi(1),
                                 domain.loVect(), domain.hiVect(), false, false);
            ylo(i,j,k,n) = lo;
            yhi(i,j,k,n) = hi;

            constexpr Real eps = 1e-6;

            Real st = (vad >= 0.) ? lo : hi;
            Real fu = (std::abs(vad) < eps) ? 0.0 : 1.0;
            Imy(i, j, k, n) = fu*st + (1.0 - fu)*0.5*(hi + lo); // store yedge
        },
        zebox, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            Real lo, hi;
            if (l_use_forces_in_trans) {
                // xxxxx lo = Ipz(i,j,k-1,n) + 0.5*l_dt*tf(i,j,k-1,n);
                // hi = Imz(i,j,k,n)   + 0.5*l_dt*tf(i,j,k-1,n);
            } else {
                lo = Ipz(i,j,k-1,n);
                hi = Imz(i,j,k  ,n);
            }

            Real wad = w_ad(i,j,k);
            auto bc = pbc[n];
            Godunov_trans_zbc_lo(i, j, k, n, q, lo, hi, wad, bc.lo(2),
                                 domain.loVect(), domain.hiVect(), false, false);
            Godunov_trans_zbc_hi(i, j, k, n, q, lo, hi, wad, bc.hi(2),
                                 domain.loVect(), domain.hiVect(), false, false);

            zlo(i,j,k,n) = lo;
            zhi(i,j,k,n) = hi;

            constexpr Real eps = 1e-6;

            Real st = (wad >= 0.) ? lo : hi;
            Real fu = (std::abs(wad) < eps) ? 0.0 : 1.0;
            Imz(i, j, k, n) = fu*st + (1.0 - fu)*0.5*(hi + lo); // store zedge
        });

    Array4<Real> xedge = Imx;
    Array4<Real> yedge = Imy;
    Array4<Real> zedge = Imz;

    Array4<Real> divu = makeArray4(Ipx.dataPtr(), grow(bx,1), 1);
    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        divu(i,j,k) = 0.0;
    });

    // We can reuse the space in Ipy and Ipz.

    //
    // X-Flux
    //
    Box xbxtmp = Box(xbx).enclosedCells().grow(0,1);
    Array4<Real> yzlo = makeArray4(Ipy.dataPtr(), amrex::surroundingNodes(xbxtmp,1),ncomp);
    Array4<Real> zylo = makeArray4(Ipz.dataPtr(), amrex::surroundingNodes(xbxtmp,2),ncomp);
    // Add d/dy term to z-faces
    // Start with {zlo,zhi} --> {zylo, zyhi} and upwind using w_ad to {zylo}
    // Add d/dz to y-faces
    // Start with {ylo,yhi} --> {yzlo, yzhi} and upwind using v_ad to {yzlo}
    amrex::ParallelFor(
    Box(zylo), ncomp,
    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        const auto bc = pbc[n]; 
        Real l_zylo, l_zyhi;
        Godunov_corner_couple_zy(l_zylo, l_zyhi,
                                 i, j, k, n, l_dt, dy, iconserv[n],
                                 zlo(i,j,k,n), zhi(i,j,k,n),
                                 q, divu, v_ad, yedge);

        Real wad = w_ad(i,j,k);
        Godunov_trans_zbc_lo(i, j, k, n, q, l_zylo, l_zyhi, wad, bc.lo(2),
                             domain.loVect(), domain.hiVect(), false, true);  
        Godunov_trans_zbc_hi(i, j, k, n, q, l_zylo, l_zyhi, wad, bc.hi(2),
                             domain.loVect(), domain.hiVect(), false, true);  

        constexpr Real eps = 1e-6;

        Real st = (wad >= 0.) ? l_zylo : l_zyhi;
        Real fu = (std::abs(wad) < eps) ? 0.0 : 1.0;
        zylo(i,j,k,n) = fu*st + (1.0 - fu) *0.5 * (l_zyhi + l_zylo);
    },
    Box(yzlo), ncomp,
    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        const auto bc = pbc[n];
        Real l_yzlo, l_yzhi;
        Godunov_corner_couple_yz(l_yzlo, l_yzhi,
                                 i, j, k, n, l_dt, dz, iconserv[n],
                                 ylo(i,j,k,n), yhi(i,j,k,n),
                                 q, divu, w_ad, zedge);

        Real vad = v_ad(i,j,k);
        Godunov_trans_ybc_lo(i, j, k, n, q, l_yzlo, l_yzhi, vad, bc.lo(1),
                             domain.loVect(), domain.hiVect(), false, true);
        Godunov_trans_ybc_hi(i, j, k, n, q, l_yzlo, l_yzhi, vad, bc.hi(1),
                             domain.loVect(), domain.hiVect(), false, true);

        constexpr Real eps = 1e-6;

        Real st = (vad >= 0.) ? l_yzlo : l_yzhi;
        Real fu = (std::abs(vad) < eps) ? 0.0 : 1.0;
        yzlo(i,j,k,n) = fu*st + (1.0 - fu) *0.5 * (l_yzhi + l_yzlo);
    });
    //
    amrex::ParallelFor(xbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        constexpr int n = 0;
        auto bc = pbc[n];
        Real stl = xlo(i,j,k,n) - (0.25*l_dt/dy)*(v_ad(i-1,j+1,k)+v_ad(i-1,j,k))*
                                  (yzlo(i-1,j+1,k,n) - yzlo(i-1,j,k,n))
                                - (0.25*l_dt/dz)*(w_ad(i-1,j,k+1)+w_ad(i-1,j,k))*
                                  (zylo(i-1,j,k+1,n) - zylo(i-1,j,k,n));
        Real sth = xhi(i,j,k,n) - (0.25*l_dt/dy)*(v_ad(i,j+1,k)+v_ad(i,j,k))*
                                  (yzlo(i,j+1,k,n) - yzlo(i,j,k,n))
                                - (0.25*l_dt/dz)*(w_ad(i,j,k+1)+w_ad(i,j,k))*
                                  (zylo(i,j,k+1,n) - zylo(i,j,k,n));

        if (!l_use_forces_in_trans)
        {
            // xxxxx stl += 0.5 * l_dt * tf(i-1,j,k,n);
            // sth += 0.5 * l_dt * tf(i  ,j,k,n);
        }

        constexpr Real eps = 1e-6;

        Godunov_cc_xbc(i, j, k, n, q, stl, sth, u_ad, bc.lo(0), bc.hi(0), dlo.x, dhi.x);

        Real st = ( (stl+sth) >= 0.) ? stl : sth;
        bool ltm = ( (stl <= 0. && sth >= 0.) || (std::abs(stl+sth) < eps) );
        qx(i,j,k) = ltm ? 0. : st;
    });
}

void incflo::predict_godunov (int lev, Real time, MultiFab& u_mac, MultiFab& v_mac,
                              MultiFab& w_mac, MultiFab const& vel, MultiFab const& rho)
{
    const int ncomp = AMREX_SPACEDIM;
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    {
        FArrayBox scratch;
        for (MFIter mfi(vel,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            Box const& bx = mfi.tilebox();
            Box const& bxg1 = amrex::grow(bx,1);
            Box const& xbx = mfi.nodaltilebox(0);
            Box const& ybx = mfi.nodaltilebox(1);
            Box const& zbx = mfi.nodaltilebox(2);

            Array4<Real> const& a_umac = u_mac.array(mfi);
            Array4<Real> const& a_vmac = v_mac.array(mfi);
            Array4<Real> const& a_wmac = w_mac.array(mfi);
            Array4<Real const> const& a_vel = vel.const_array(mfi);

            scratch.resize(bxg1, ncomp*12+3);
//            Elixir eli = scratch.elixir(); // not needed because of streamSynchronize later
            Real* p = scratch.dataPtr();

            Array4<Real> Imx = makeArray4(p,bxg1,ncomp);
            p +=         Imx.size();
            Array4<Real> Ipx = makeArray4(p,bxg1,ncomp);
            p +=         Ipx.size();
            Array4<Real> Imy = makeArray4(p,bxg1,ncomp);
            p +=         Imy.size();
            Array4<Real> Ipy = makeArray4(p,bxg1,ncomp);
            p +=         Ipy.size();
            Array4<Real> Imz = makeArray4(p,bxg1,ncomp);
            p +=         Imz.size();
            Array4<Real> Ipz = makeArray4(p,bxg1,ncomp);
            p +=         Ipz.size();
            Array4<Real> u_ad = makeArray4(p,Box(bx).grow(1,1).grow(2,1).surroundingNodes(0),1);
            p +=         u_ad.size();
            Array4<Real> v_ad = makeArray4(p,Box(bx).grow(0,1).grow(2,1).surroundingNodes(1),1);
            p +=         v_ad.size();
            Array4<Real> w_ad = makeArray4(p,Box(bx).grow(0,1).grow(1,1).surroundingNodes(2),1);
            p +=         w_ad.size();

            make_ppm_integrals (lev, bxg1, AMREX_SPACEDIM, Imx, Ipx, Imy, Ipy, Imz, Ipz, a_vel, a_vel);

            make_trans_velocities(lev, Box(u_ad), Box(v_ad), Box(w_ad),
                                  u_ad, v_ad, w_ad,
                                  Imx, Ipx, Imy, Ipy, Imz, Ipz, a_vel);

            predict_godunov_on_box(lev, bx, ncomp, xbx, ybx, zbx, a_umac, a_vmac, a_wmac,
                                   a_vel, u_ad, v_ad, w_ad, Imx, Ipx, Imy, Ipy, Imz, Ipz, p);

            Gpu::streamSynchronize();  // otherwise we might be using too much memory
        }
    }
     
    VisMF::Write(u_mac, "u_mac");

    amrex::Abort("end of predict_godunov");
}

void
incflo::incflo_predict_godunov ( int lev, Real time,
                                 Vector< std::unique_ptr<MultiFab> >& vel_in,
                                 Vector< std::unique_ptr<MultiFab> >& vel_forces_in)
{
    BL_PROFILE("incflo::incflo_predict_godunov");

    Box domain(geom[lev].Domain());

    BCRec dom_bc;
    {
      // const int* lo_bc = phys_bc.lo();
      // const int* hi_bc = phys_bc.hi();
      // HACK -- just set to all int_dir as stand-in for periodic
      dom_bc.setLo(0,BCType::int_dir);
      dom_bc.setHi(0,BCType::int_dir);
      dom_bc.setLo(1,BCType::int_dir);
      dom_bc.setHi(1,BCType::int_dir);
      dom_bc.setLo(2,BCType::int_dir);
      dom_bc.setHi(2,BCType::int_dir);
    }

    for (MFIter mfi(*vel_in[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        // Tilebox
        Box bx = mfi.tilebox();

        Gpu::ManagedVector<BCRec> bc(3);
        for (int n = 0; n < 3; ++n)
            setBC(bx, geom[lev].Domain(), dom_bc, bc[n]);

        // Create temporaries to hold the normal velocities on faces that we will *only* use for upwinding
        // in the computation of the transverse derivatives
        FArrayBox uad((*m_u_mac[lev])[mfi].box(), 1);
        Elixir uad_eli = uad.elixir();
        FArrayBox vad((*m_v_mac[lev])[mfi].box(), 1);
        Elixir vad_eli = vad.elixir();
        FArrayBox wad((*m_w_mac[lev])[mfi].box(), 1);
        Elixir wad_eli = wad.elixir();

        // The output of this call is predicted normal velocities on faces that we will *only* use for upwinding
        // in the computation of the transverse derivatives
        incflo_make_trans_velocities(lev, bx, (*vel_in[lev]).array(mfi), 
                                     uad.array(), vad.array(), wad.array(), 
                                     (*vel_forces_in[lev]).array(mfi), bc);

        // The output of this call is the predicted normal velocities on faces -- we store these in 
        // {m_u_mac, m_v_mac, m_w_mac} but they have not yet been projected
        incflo_predict_godunov_on_box(lev, bx, (*vel_in[lev]).array(mfi), 
                                      uad.array(), vad.array(), wad.array(), 
                                      m_u_mac[lev]->array(mfi), m_v_mac[lev]->array(mfi), m_w_mac[lev]->array(mfi),
                                      (*vel_forces_in[lev]).array(mfi), bc);
    }
}

void
incflo::incflo_predict_godunov_on_box (const int lev, Box& bx,
                                       const Array4<Real> &a_vel,
                                       const Array4<const Real> &u_ad, 
                                       const Array4<const Real> &v_ad, 
                                       const Array4<const Real> &w_ad,
                                       const Array4<      Real> &u_face, 
                                       const Array4<      Real> &v_face, 
                                       const Array4<      Real> &w_face,
                                       const Array4<Real> &tf, 
                                       const Gpu::ManagedVector<BCRec> &BCs)
{
    Box domain(geom[lev].Domain());

    // Need to do these to make GPUs happy since dt and use_forces_in_trans are members of incflo 
    Real l_dt = dt;
    Real l_use_forces_in_trans = use_forces_in_trans;

    // We only predict the normal velocity in this routine but we need transverse terms
    //   which means we need all velocities on all faces
    int ncomp = 3;

    // We only predict velocity in this routine, which is convectively not conservatively, updated
    GpuArray<int,3> iconserv{0,0,0};

    auto const g2bx = amrex::grow(bx, 2); 
    auto const g1bx = amrex::grow(bx, 1); 

    const Real dx = geom[lev].CellSize(0); 
    const Real dy = geom[lev].CellSize(1); 
    const Real dz = geom[lev].CellSize(2); 
   
    auto const gxbx = amrex::grow(g2bx,0, 1); 
    FArrayBox Imxf(gxbx, ncomp); 
    FArrayBox Ipxf(gxbx, ncomp);
    Elixir Imxeli = Imxf.elixir(); 
    Elixir Ipxeli = Ipxf.elixir(); 
    auto const Imx = Imxf.array(); 
    auto const Ipx = Ipxf.array(); 

    auto const gybx = amrex::grow(g2bx,1, 1); 
    FArrayBox Imyf(gybx, ncomp); 
    FArrayBox Ipyf(gybx, ncomp);
    Elixir Imyeli = Imyf.elixir(); 
    Elixir Ipyeli = Ipyf.elixir(); 
    auto const Imy = Imyf.array(); 
    auto const Ipy = Ipyf.array(); 
   
    auto const gzbx = amrex::grow(g2bx,2, 1); 
    FArrayBox Imzf(gzbx, ncomp); 
    FArrayBox Ipzf(gzbx, ncomp);
    Elixir Imzeli = Imzf.elixir(); 
    Elixir Ipzeli = Ipzf.elixir(); 
    auto const Imz = Imzf.array(); 
    auto const Ipz = Ipzf.array(); 

    /* Temporary Edge States */ 
    FArrayBox xedgef(g2bx, ncomp); 
    FArrayBox yedgef(g2bx, ncomp); 
    FArrayBox zedgef(g2bx, ncomp); 

    Elixir xedeli = xedgef.elixir(); 
    Elixir yedeli = yedgef.elixir(); 
    Elixir zedeli = zedgef.elixir(); 

    auto const xedge = xedgef.array(); 
    auto const yedge = yedgef.array(); 
    auto const zedge = zedgef.array(); 

    auto const xgbx = surroundingNodes(g1bx, 0);
    auto const ygbx = surroundingNodes(g1bx, 1); 
    auto const zgbx = surroundingNodes(g1bx, 2); 

    BCRec const* pbc = BCs.data();

    AMREX_PARALLEL_FOR_4D (g1bx, ncomp, i, j, k, n, {
        const auto bc = pbc[n];
        Godunov_ppm_pred_x(i, j, k, n, l_dt, dx, a_vel, a_vel, Imx, Ipx, bc,
                         domain.loVect()[0], domain.hiVect()[0]);
    });

    AMREX_PARALLEL_FOR_4D (g1bx, ncomp, i, j, k, n, {
        const auto bc = pbc[n];
        Godunov_ppm_pred_y(i, j, k, n, l_dt, dy, a_vel, a_vel, Imy, Ipy, bc,
                         domain.loVect()[1], domain.hiVect()[1]);
    });

    AMREX_PARALLEL_FOR_4D (g1bx, ncomp, i, j, k, n, {
        const auto bc = pbc[n];
        Godunov_ppm_pred_z(i, j, k, n, l_dt, dz, a_vel, a_vel, Imz, Ipz, bc,
                         domain.loVect()[2], domain.hiVect()[2]);
    }); 

    FArrayBox xlf(xgbx, ncomp); 
    FArrayBox xhf(xgbx, ncomp);
    Elixir xleli = xlf.elixir(); 
    Elixir xheli = xhf.elixir(); 
    auto const xlo = xlf.array(); 
    auto const xhi = xhf.array(); 

    FArrayBox ylf(ygbx, ncomp); 
    FArrayBox yhf(ygbx, ncomp);
    Elixir yleli = ylf.elixir(); 
    Elixir yheli = yhf.elixir(); 
    auto const ylo = ylf.array(); 
    auto const yhi = yhf.array(); 
   
    FArrayBox zlf(zgbx, ncomp); 
    FArrayBox zhf(zgbx, ncomp);
    Elixir zleli = zlf.elixir(); 
    Elixir zheli = zhf.elixir(); 
    auto const zlo = zlf.array(); 
    auto const zhi = zhf.array();

    auto const txbx = surroundingNodes(grow(g1bx, 0, -1), 0);
    auto const tybx = surroundingNodes(grow(g1bx, 1, -1), 1); 
    auto const tzbx = surroundingNodes(grow(g1bx, 2, -1), 2);
    
// --------------------X -------------------------------------------------
    AMREX_PARALLEL_FOR_4D(txbx, ncomp, i, j, k, n, 
    {
        Real lo;
        Real hi;

        if (l_use_forces_in_trans) {
           lo = Ipx(i-1,j,k,n) + 0.5*l_dt*tf(i-1,j,k,n);
           hi = Imx(i  ,j,k,n) + 0.5*l_dt*tf(i  ,j,k,n);
        } else {
           lo = Ipx(i-1,j,k,n);
           hi = Imx(i  ,j,k,n);
        }

        Real uad = u_ad(i,j,k);
        auto  bc = pbc[n];  

        Godunov_trans_xbc_lo(i, j, k, n, a_vel, lo, hi, uad, bc.lo(0),
                             domain.loVect(), domain.hiVect(), false, false);  
        Godunov_trans_xbc_hi(i, j, k, n, a_vel, lo, hi, uad, bc.hi(0),
                             domain.loVect(), domain.hiVect(), false, false);  
        xlo(i,j,k,n) = lo; 
        xhi(i,j,k,n) = hi;

        Real eps = 1e-6;

        Real st = (uad >= 0.) ? lo : hi;
        Real fu = (std::abs(uad) < eps) ? 0.0 : 1.0; 
        xedge(i, j, k, n) = fu*st + (1.0 - fu) *0.5 * (hi + lo);
    }); 

    Imxeli.clear(); 
    Ipxeli.clear(); 

//--------------------Y -------------------------------------------------
     AMREX_PARALLEL_FOR_4D(tybx, ncomp, i, j, k, n, 
     {
        Real lo;
        Real hi;

        if (l_use_forces_in_trans) {
           lo = Ipy(i,j-1,k,n) + 0.5*l_dt*tf(i,j-1,k,n);
           hi = Imy(i,j,k,n)   + 0.5*l_dt*tf(i,j-1,k,n);
        } else {
           lo = Ipy(i,j-1,k,n);
           hi = Imy(i,j  ,k,n);
        }

        Real vad = v_ad(i,j,k);
        auto  bc = pbc[n];  

        Godunov_trans_ybc_lo(i, j, k, n, a_vel, lo, hi, vad, bc.lo(1), 
                             domain.loVect(), domain.hiVect(), false, false);  
        Godunov_trans_ybc_hi(i, j, k, n, a_vel, lo, hi, vad, bc.hi(1),
                             domain.loVect(), domain.hiVect(), false, false);  
        ylo(i,j,k,n) = lo; 
        yhi(i,j,k,n) = hi; 

        Real eps = 1e-6;

        Real st = (vad >= 0.) ? lo : hi;
        Real fu = (std::abs(vad) < eps) ? 0.0 : 1.0; 
        yedge(i, j, k, n) = fu*st + (1.0 - fu)*0.5*(hi + lo);
    });

    Imyeli.clear(); 
    Ipyeli.clear(); 

//------------------ Z ---------------------------------------------------
     AMREX_PARALLEL_FOR_4D(tzbx, ncomp, i, j, k, n, 
     {
        Real lo;
        Real hi;

        if (l_use_forces_in_trans) {
           lo = Ipz(i,j,k-1,n) + 0.5*l_dt*tf(i,j,k-1,n);
           hi = Imz(i,j,k,n)   + 0.5*l_dt*tf(i,j,k-1,n);
        } else {
           lo = Ipz(i,j,k-1,n);
           hi = Imz(i,j,k  ,n);
        }

        Real wad = w_ad(i,j,k);
        auto  bc = pbc[n];  
        Godunov_trans_zbc_lo(i, j, k, n, a_vel, lo, hi, wad, bc.lo(2),
                             domain.loVect(), domain.hiVect(), false, false);  
        Godunov_trans_zbc_hi(i, j, k, n, a_vel, lo, hi, wad, bc.hi(2),
                             domain.loVect(), domain.hiVect(), false, false);  

        zlo(i,j,k,n) = lo; 
        zhi(i,j,k,n) = hi; 

        Real eps = 1e-6;
        
        Real st = (wad >= 0.) ? lo : hi;
        Real fu = (std::abs(wad) < eps) ? 0.0 : 1.0; 

        zedge(i, j, k, n) = fu*st + (1.0 - fu)*0.5*(hi + lo);
    });     

    Imzeli.clear(); 
    Ipzeli.clear(); 
 

//-----------------Create temporary fabs for corner/transverse 

    //X fabs
    FArrayBox xylof(g1bx, ncomp); 
    FArrayBox xyhif(g1bx, ncomp); 
    FArrayBox xzlof(g1bx, ncomp); 
    FArrayBox xzhif(g1bx, ncomp); 
    Elixir xyle = xylof.elixir(); 
    Elixir xyhe = xyhif.elixir(); 
    Elixir xzle = xzlof.elixir(); 
    Elixir xzhe = xzhif.elixir(); 
    const auto xylo = xylof.array(); 
    const auto xyhi = xyhif.array(); 
    const auto xzlo = xzlof.array(); 
    const auto xzhi = xzhif.array(); 

    //Y fabs
    FArrayBox yxlof(g1bx, ncomp); 
    FArrayBox yxhif(g1bx, ncomp); 
    FArrayBox yzlof(g1bx, ncomp); 
    FArrayBox yzhif(g1bx, ncomp); 
    Elixir yxle = yxlof.elixir(); 
    Elixir yxhe = yxhif.elixir(); 
    Elixir yzle = yzlof.elixir(); 
    Elixir yzhe = yzhif.elixir(); 
    const auto yxlo = yxlof.array(); 
    const auto yxhi = yxhif.array(); 
    const auto yzlo = yzlof.array(); 
    const auto yzhi = yzhif.array(); 

    //Z fabs 
    FArrayBox zxlof(g1bx, ncomp); 
    FArrayBox zxhif(g1bx, ncomp); 
    FArrayBox zylof(g1bx, ncomp); 
    FArrayBox zyhif(g1bx, ncomp); 
    Elixir zxle = zxlof.elixir(); 
    Elixir zxhe = zxhif.elixir(); 
    Elixir zyle = zylof.elixir(); 
    Elixir zyhe = zyhif.elixir(); 
    const auto zxlo = zxlof.array(); 
    const auto zxhi = zxhif.array(); 
    const auto zylo = zylof.array(); 
    const auto zyhi = zyhif.array(); 

    auto const xybx = surroundingNodes(grow(bx, 2, 1), 0);
    auto const xzbx = surroundingNodes(grow(bx, 1, 1), 0);  
    auto const yxbx = surroundingNodes(grow(bx, 2, 1), 1); 
    auto const yzbx = surroundingNodes(grow(bx, 0, 1), 1); 
    auto const zxbx = surroundingNodes(grow(bx, 1, 1), 2); 
    auto const zybx = surroundingNodes(grow(bx, 0, 1), 2); 

    // We need divu_fab only to pass to Godunov_corner_couple but it isn't actually used
    FArrayBox divu_fab(grow(bx,2), 1);
    divu_fab.setVal(0.0);
    Elixir divu_eli = divu_fab.elixir();
    const auto divu_arr = divu_fab.array(); 

/*------------------Now perform corner coupling */
    // Add d/dx term to y-faces and upwind using v_ad to {yxlo}
    // Start with {ylo,yhi} --> {yxlo, yxhi}
    AMREX_PARALLEL_FOR_4D (yxbx, ncomp, i, j, k, n, 
    {
        const auto bc = pbc[n];
        Godunov_corner_couple_yx(yxlo(i,j,k,n), yxhi(i,j,k,n),
                                 i, j, k, n, l_dt, dx, iconserv[n],
                                 ylo(i,j,k,n), yhi(i,j,k,n),
                                 a_vel, divu_arr, u_ad, xedge);

        Real vad = v_ad(i,j,k);
        Godunov_trans_ybc_lo(i, j, k, n, a_vel, yxlo(i,j,k,n), yxhi(i,j,k,n), vad, bc.lo(1),
                             domain.loVect(), domain.hiVect(), true, false);  
        Godunov_trans_ybc_hi(i, j, k, n, a_vel, yxlo(i,j,k,n), yxhi(i,j,k,n), vad, bc.hi(1),
                             domain.loVect(), domain.hiVect(), true, false);  

        Real eps = 1e-6;

        Real st = (vad >= 0.) ? yxlo(i,j,k,n) : yxhi(i,j,k,n);
        Real fu = (std::abs(vad) < eps) ? 0.0 : 1.0; 
        yxlo(i,j,k,n) = fu*st + (1.0 - fu) *0.5 * (yxhi(i,j,k,n) + yxlo(i,j,k,n));
    }); 
    
    // Add d/dx term to z-faces
    // Start with {zlo,zhi} --> {zxlo, zxhi} and upwind using w_ad to {zxlo}
    AMREX_PARALLEL_FOR_4D (zxbx, ncomp, i, j, k, n, 
    {
        const auto bc = pbc[n]; 
        Godunov_corner_couple_zx(zxlo(i,j,k,n), zxhi(i,j,k,n),
                                 i, j, k, n, l_dt, dx, iconserv[n],
                                 zlo(i,j,k,n), zhi(i,j,k,n),
                                 a_vel, divu_arr, u_ad, xedge);

        Real wad = w_ad(i,j,k);
        Godunov_trans_zbc_lo(i, j, k, n, a_vel, zxlo(i,j,k,n), zxhi(i,j,k,n), wad, bc.lo(2), 
                             domain.loVect(), domain.hiVect(), true, false);  
        Godunov_trans_zbc_hi(i, j, k, n, a_vel, zxlo(i,j,k,n), zxhi(i,j,k,n), wad, bc.hi(2), 
                             domain.loVect(), domain.hiVect(), true, false);  

        Real eps = 1e-6;

        Real st = (wad >= 0.) ? zxlo(i,j,k,n) : zxhi(i,j,k,n);
        Real fu = (std::abs(wad) < eps) ? 0.0 : 1.0; 
        zxlo(i,j,k,n) = fu*st + (1.0 - fu) *0.5 * (zxhi(i,j,k,n) + zxlo(i,j,k,n));
    });

    // Add d/dy term to x-faces
    // Start with {xlo,xhi} --> {xylo, xyhi} and upwind using u_ad to {xylo}
    AMREX_PARALLEL_FOR_4D (xybx, ncomp, i, j, k, n, 
    {
        const auto bc = pbc[n]; 
        Godunov_corner_couple_xy(xylo(i,j,k,n), xyhi(i,j,k,n),
                                 i, j, k, n, l_dt, dy, iconserv[n],
                                 xlo(i,j,k,n), xhi(i,j,k,n),
                                 a_vel, divu_arr, v_ad, yedge);

        Real uad = u_ad(i,j,k);
        Godunov_trans_xbc_lo(i, j, k, n, a_vel, xylo(i,j,k,n), xyhi(i,j,k,n), uad, bc.lo(0),
                             domain.loVect(), domain.hiVect(), true, false);  
        Godunov_trans_xbc_hi(i, j, k, n, a_vel, xylo(i,j,k,n), xyhi(i,j,k,n), uad, bc.hi(0),
                             domain.loVect(), domain.hiVect(), true, false);  

        Real eps = 1e-6;

        Real st = (uad >= 0.) ? xylo(i,j,k,n) : xyhi(i,j,k,n);
        Real fu = (std::abs(uad) < eps) ? 0.0 : 1.0; 
        xylo(i,j,k,n) = fu*st + (1.0 - fu) *0.5 * (xyhi(i,j,k,n) + xylo(i,j,k,n));

    }); 

    // Add d/dy term to z-faces
    // Start with {zlo,zhi} --> {zylo, zyhi} and upwind using w_ad to {zylo}
    AMREX_PARALLEL_FOR_4D (zybx, ncomp, i, j, k, n, 
    {
        const auto bc = pbc[n]; 
        Godunov_corner_couple_zy(zylo(i,j,k,n), zyhi(i,j,k,n),
                                 i, j, k, n, l_dt, dy, iconserv[n],
                                 zlo(i,j,k,n), zhi(i,j,k,n),
                                 a_vel, divu_arr, v_ad, yedge);

        Real wad = w_ad(i,j,k);
        Godunov_trans_zbc_lo(i, j, k, n, a_vel, zylo(i,j,k,n), zyhi(i,j,k,n), wad, bc.lo(2),
                             domain.loVect(), domain.hiVect(), false, true);  
        Godunov_trans_zbc_hi(i, j, k, n, a_vel, zylo(i,j,k,n), zyhi(i,j,k,n), wad, bc.hi(2),
                             domain.loVect(), domain.hiVect(), false, true);  

        Real eps = 1e-6;

        Real st = (wad >= 0.) ? zylo(i,j,k,n) : zyhi(i,j,k,n);
        Real fu = (std::abs(wad) < eps) ? 0.0 : 1.0; 
        zylo(i,j,k,n) = fu*st + (1.0 - fu) *0.5 * (zyhi(i,j,k,n) + zylo(i,j,k,n));

    });

    // Add d/dz to x-faces
    // Start with {xlo,xhi} --> {xzlo, xzhi} and upwind using u_ad to {xzlo}
    AMREX_PARALLEL_FOR_4D (xzbx, ncomp, i, j, k, n, 
    {
        const auto bc = pbc[n]; 
        Godunov_corner_couple_xz(xzlo(i,j,k,n), xzhi(i,j,k,n),
                                 i, j, k, n, l_dt, dz, iconserv[n],
                                 xlo(i,j,k,n), xhi(i,j,k,n),
                                 a_vel, divu_arr, w_ad, zedge);

        Real uad = u_ad(i,j,k);
        Godunov_trans_xbc_lo(i, j, k, n, a_vel, xzlo(i,j,k,n), xzhi(i,j,k,n), uad, bc.lo(0),
                             domain.loVect(), domain.hiVect(), false, true);  
        Godunov_trans_xbc_hi(i, j, k, n, a_vel, xzlo(i,j,k,n), xzhi(i,j,k,n), uad, bc.hi(0),
                             domain.loVect(), domain.hiVect(), false, true);  

        Real eps = 1e-6;

        Real st = (uad >= 0.) ? xzlo(i,j,k,n) : xzhi(i,j,k,n);
        Real fu = (std::abs(uad) < eps) ? 0.0 : 1.0; 
        xzlo(i,j,k,n) = fu*st + (1.0 - fu) *0.5 * (xzhi(i,j,k,n) + xzlo(i,j,k,n));
    });  

    // Add d/dz to y-faces
    // Start with {ylo,yhi} --> {yzlo, yzhi} and upwind using v_ad to {yzlo}
    AMREX_PARALLEL_FOR_4D (yzbx, ncomp, i, j, k, n, 
    {
        const auto bc = pbc[n]; 
        Godunov_corner_couple_yz(yzlo(i,j,k,n), yzhi(i,j,k,n),
                                 i, j, k, n, l_dt, dz, iconserv[n],
                                 ylo(i,j,k,n), yhi(i,j,k,n),
                                 a_vel, divu_arr, w_ad, zedge);

        Real vad = v_ad(i,j,k);
        Godunov_trans_ybc_lo(i, j, k, n, a_vel, yzlo(i,j,k,n), yzhi(i,j,k,n), vad, bc.lo(1),
                             domain.loVect(), domain.hiVect(), false, true);  
        Godunov_trans_ybc_hi(i, j, k, n, a_vel, yzlo(i,j,k,n), yzhi(i,j,k,n), vad, bc.hi(1),
                             domain.loVect(), domain.hiVect(), false, true);  

        Real eps = 1e-6;

        Real st = (vad >= 0.) ? yzlo(i,j,k,n) : yzhi(i,j,k,n);
        Real fu = (std::abs(vad) < eps) ? 0.0 : 1.0; 
        yzlo(i,j,k,n) = fu*st + (1.0 - fu) *0.5 * (yzhi(i,j,k,n) + yzlo(i,j,k,n));
    });

   const auto xbx = surroundingNodes(bx, 0); 
   const auto ybx = surroundingNodes(bx, 1); 
   const auto zbx = surroundingNodes(bx, 2);
 
//--------------------------------------- X -------------------------------------- 
    /* Final Update of Faces */ 
    AMREX_PARALLEL_FOR_3D(xbx, i, j, k,
    {
        // We only care about x-velocity on x-faces here
        int n = 0;

        auto bc = pbc[n]; 
        Real stl = xlo(i,j,k,n) - (0.25*l_dt/dy)*(v_ad(i-1,j+1,k)+v_ad(i-1,j,k))*
                                  (yzlo(i-1,j+1,k,n) - yzlo(i-1,j,k,n))
                                - (0.25*l_dt/dz)*(w_ad(i-1,j,k+1)+w_ad(i-1,j,k))*
                                  (zylo(i-1,j,k+1,n) - zylo(i-1,j,k,n));

        Real sth = xhi(i,j,k,n) - (0.25*l_dt/dy)*(v_ad(i,j+1,k)+v_ad(i,j,k))*
                                  (yzlo(i,j+1,k,n) - yzlo(i,j,k,n))
                                - (0.25*l_dt/dz)*(w_ad(i,j,k+1)+w_ad(i,j,k))*
                                  (zylo(i,j,k+1,n) - zylo(i,j,k,n));

        if (!l_use_forces_in_trans)
        { 
           stl += 0.5 * l_dt * tf(i-1,j,k,n);
           sth += 0.5 * l_dt * tf(i  ,j,k,n);
        }

        Real eps = 1e-6;
       
        Godunov_cc_xbc(i, j, k, n, a_vel, stl, sth, u_ad, bc.lo(0), bc.hi(0),
                       domain.loVect()[0], domain.hiVect()[0]);  

        Real st = ( (stl+sth) >= 0.) ? stl : sth;
        bool ltm = ( (stl <= 0. && sth >= 0.) || (std::abs(stl+sth) < eps) );
        u_face(i,j,k) = ltm ? 0. : st;
    }); 

//-------------------------------------- Y ------------------------------------            
    AMREX_PARALLEL_FOR_3D(ybx, i, j, k,
    {
        // We only care about y-velocity on y-faces here
        int n = 1;

        auto bc = pbc[n]; 

        Real stl = ylo(i,j,k,n) - (0.25*l_dt/dx)*(u_ad(i+1,j-1,k)+u_ad(i,j-1,k))*
                                  (xzlo(i+1,j-1,k,n) - xzlo(i,j-1,k,n))
                                - (0.25*l_dt/dz)*(w_ad(i,j-1,k+1)+w_ad(i,j-1,k))*
                                  (zxlo(i,j-1,k+1,n) - zxlo(i,j-1,k,n));

        Real sth = yhi(i,j,k,n) - (0.25*l_dt/dx)*(u_ad(i+1,j,k)+u_ad(i,j,k))*
                                  (xzlo(i+1,j,k,n) - xzlo(i,j,k,n))
                                - (0.25*l_dt/dz)*(w_ad(i,j,k+1)+w_ad(i,j,k))*
                                  (zxlo(i,j,k+1,n) - zxlo(i,j,k,n));

        if (!l_use_forces_in_trans)
        { 
           stl += 0.5 * l_dt * tf(i,j-1,k,n);
           sth += 0.5 * l_dt * tf(i,j  ,k,n);
        }

        Real eps = 1e-6;

        Godunov_cc_ybc(i, j, k, n, a_vel, stl, sth, v_ad, bc.lo(1), bc.hi(1), 
                       domain.loVect()[1], domain.hiVect()[1]); 

        Real st = ( (stl+sth) >= 0.) ? stl : sth;
        bool ltm = ( (stl <= 0. && sth >= 0.) || (std::abs(stl+sth) < eps) );
        v_face(i,j,k) = ltm ? 0. : st;
     });

//----------------------------------- Z ----------------------------------------- 
     AMREX_PARALLEL_FOR_3D(zbx, i, j, k,
     {
        // We only care about z-velocity on z-faces here
        int n = 2;

        auto bc = pbc[n]; 

        Real stl = zlo(i,j,k,n) - (0.25*l_dt/dx)*(u_ad(i+1,j,k-1)+u_ad(i,j,k-1))*
                                  (xylo(i+1,j,k-1,n) - xylo(i,j,k-1,n))
                                - (0.25*l_dt/dy)*(v_ad(i,j+1,k-1)+v_ad(i,j,k-1))*
                                  (yxlo(i,j+1,k-1,n) - yxlo(i,j,k-1,n));
     
        Real sth = zhi(i,j,k,n) - (0.25*l_dt/dx)*(u_ad(i+1,j,k)+u_ad(i,j,k))*
                                  (xylo(i+1,j,k,n) - xylo(i,j,k,n))
                                - (0.25*l_dt/dy)*(v_ad(i,j+1,k)+v_ad(i,j,k))*
                                  (yxlo(i,j+1,k,n) - yxlo(i,j,k,n));

        if (!l_use_forces_in_trans)
        { 
           stl += 0.5 * l_dt * tf(i,j,k-1,n);
           sth += 0.5 * l_dt * tf(i,j,k  ,n);
        }

        Godunov_cc_zbc(i, j, k, n, a_vel, stl, sth, w_ad, bc.lo(2), bc.hi(2), 
                       domain.loVect()[2], domain.hiVect()[2]);  

        Real eps = 1e-6;

        Real st = ( (stl+sth) >= 0.) ? stl : sth;
        bool ltm = ( (stl <= 0. && sth >= 0.) || (std::abs(stl+sth) < eps) );
        w_face(i,j,k) = ltm ? 0. : st;
    }); 
}

void
incflo::incflo_make_trans_velocities (const int lev, Box& bx,
                                      const Array4<const Real> &a_vel,
                                      const Array4<      Real> &u_ad, 
                                      const Array4<      Real> &v_ad, 
                                      const Array4<      Real> &w_ad,
                                      const Array4<Real> &tf, 
                                      const Gpu::ManagedVector<BCRec> &BCs)
{
    Box domain(geom[lev].Domain());

    // Need to do these to make GPUs happy since dt and use_forces_in_trans are members of incflo 
    Real l_dt = dt;
    Real l_use_forces_in_trans = use_forces_in_trans;

    // We only predict the normal velocity in this routine but because we are calling
    //   the same Godunov_ppm_pred routine as called by ppm_predict, we must pass 
    //   arrays with all three components
    int ncomp = 3;

    auto const g2bx = amrex::grow(bx, 2); 
    auto const g1bx = amrex::grow(bx, 1); 

    const Real dx = geom[lev].CellSize(0); 
    const Real dy = geom[lev].CellSize(1); 
    const Real dz = geom[lev].CellSize(2); 
   
    auto const gxbx = amrex::grow(g2bx,0, 1); 
    FArrayBox Imxf(gxbx, ncomp); 
    FArrayBox Ipxf(gxbx, ncomp);
    Elixir Imxeli = Imxf.elixir(); 
    Elixir Ipxeli = Ipxf.elixir(); 
    auto const Imx = Imxf.array(); 
    auto const Ipx = Ipxf.array(); 

    auto const gybx = amrex::grow(g2bx,1, 1); 
    FArrayBox Imyf(gybx, ncomp); 
    FArrayBox Ipyf(gybx, ncomp);
    Elixir Imyeli = Imyf.elixir(); 
    Elixir Ipyeli = Ipyf.elixir(); 
    auto const Imy = Imyf.array(); 
    auto const Ipy = Ipyf.array(); 
   
    auto const gzbx = amrex::grow(g2bx,2, 1); 
    FArrayBox Imzf(gzbx, ncomp); 
    FArrayBox Ipzf(gzbx, ncomp);
    Elixir Imzeli = Imzf.elixir(); 
    Elixir Ipzeli = Ipzf.elixir(); 
    auto const Imz = Imzf.array(); 
    auto const Ipz = Ipzf.array(); 

    BCRec const* pbc = BCs.data();

    /* Use PPM to generate Im and Ip */
 
    AMREX_PARALLEL_FOR_4D (g1bx, ncomp, i, j, k, n, {
        const auto bc = pbc[n];
        Godunov_ppm_pred_x(i, j, k, n, l_dt, dx, a_vel, a_vel, Imx, Ipx, bc,
                         domain.loVect()[0], domain.hiVect()[0]);
    });

    AMREX_PARALLEL_FOR_4D (g1bx, ncomp, i, j, k, n, {
        const auto bc = pbc[n];
        Godunov_ppm_pred_y(i, j, k, n, l_dt, dy, a_vel, a_vel, Imy, Ipy, bc,
                         domain.loVect()[1], domain.hiVect()[1]);
    });

    AMREX_PARALLEL_FOR_4D (g1bx, ncomp, i, j, k, n, {
        const auto bc = pbc[n];
        Godunov_ppm_pred_z(i, j, k, n, l_dt, dz, a_vel, a_vel, Imz, Ipz, bc,
                         domain.loVect()[2], domain.hiVect()[2]);
    }); 

    auto const txbx = surroundingNodes(grow(g1bx, 0, -1), 0);
    auto const tybx = surroundingNodes(grow(g1bx, 1, -1), 1); 
    auto const tzbx = surroundingNodes(grow(g1bx, 2, -1), 2);

    Gpu::synchronize();
    AMREX_GPU_ERROR_CHECK();
    
// --------------------X -------------------------------------------------
    AMREX_PARALLEL_FOR_3D(txbx, i, j, k, 
    {
        // We only care about x-velocity on x-faces here
        int n = 0;

        Real lo;
        Real hi;

        if (l_use_forces_in_trans) {
           lo = Ipx(i-1,j,k,n) + 0.5*l_dt*tf(i-1,j,k,n);
           hi = Imx(i  ,j,k,n) + 0.5*l_dt*tf(i  ,j,k,n);
        } else {
           lo = Ipx(i-1,j,k,0);
           hi = Imx(i  ,j,k,0);
        }

        auto bc = pbc[n];  
        Godunov_trans_xbc_lo(i, j, k, n, a_vel, lo, hi, lo, bc.lo(0),
                             domain.loVect(), domain.hiVect(), false, false);  
        Godunov_trans_xbc_hi(i, j, k, n, a_vel, lo, hi, lo, bc.hi(0), 
                             domain.loVect(), domain.hiVect(), false, false);  
        Real eps = 1e-6;

        Real st = ( (lo+hi) >= 0.) ? lo : hi;
        bool ltm = ( (lo <= 0. && hi >= 0.) || (std::abs(lo+hi) < eps) );
        u_ad(i,j,k) = ltm ? 0. : st;
    }); 

    Imxeli.clear(); 
    Ipxeli.clear(); 

    Gpu::synchronize();
    AMREX_GPU_ERROR_CHECK();

//--------------------Y -------------------------------------------------
     AMREX_PARALLEL_FOR_3D(tybx, i, j, k,
     {
        // We only care about y-velocity on y-faces here
        int n = 1;

        Real lo;
        Real hi;

        if (l_use_forces_in_trans) {
           lo = Ipy(i,j-1,k,n) + 0.5*l_dt*tf(i,j-1,k,n);
           hi = Imy(i,j,k,n)   + 0.5*l_dt*tf(i,j-1,k,n);
        } else {
           lo = Ipy(i,j-1,k,n);
           hi = Imy(i,j  ,k,n);
        }

        auto bc = pbc[n];  
        Godunov_trans_ybc_lo(i, j, k, n, a_vel, lo, hi, lo, bc.lo(1), 
                             domain.loVect(), domain.hiVect(), false, false);  
        Godunov_trans_ybc_hi(i, j, k, n, a_vel, lo, hi, lo, bc.hi(1),
                             domain.loVect(), domain.hiVect(), false, false);  
        Real eps = 1e-6;

        Real st = ( (lo+hi) >= 0.) ? lo : hi;
        bool ltm = ( (lo <= 0. && hi >= 0.) || (std::abs(lo+hi) < eps) );
        v_ad(i,j,k) = ltm ? 0. : st;
    });

    Imyeli.clear(); 
    Ipyeli.clear(); 

//------------------ Z ---------------------------------------------------
     AMREX_PARALLEL_FOR_3D(tzbx, i, j, k,
     {
        // We only care about z-velocity on z-faces here
        int n = 2;

        Real lo;
        Real hi;

        if (l_use_forces_in_trans) {
           lo = Ipz(i,j,k-1,n) + 0.5*l_dt*tf(i,j,k-1,n);
           hi = Imz(i,j,k,n)   + 0.5*l_dt*tf(i,j,k-1,n);
        } else {
           lo = Ipz(i,j,k-1,n);
           hi = Imz(i,j,k  ,n);
        }

        auto bc = pbc[n];  
        Godunov_trans_zbc_lo(i, j, k, n, a_vel, lo, hi, lo, bc.lo(2),
                             domain.loVect(), domain.hiVect(), false, false);  
        Godunov_trans_zbc_hi(i, j, k, n, a_vel, lo, hi, lo, bc.hi(2),
                             domain.loVect(), domain.hiVect(), false, false);  

        Real eps = 1e-6;

        Real st = ( (lo+hi) >= 0.) ? lo : hi;
        bool ltm = ( (lo <= 0. && hi >= 0.) || (std::abs(lo+hi) < eps) );
        w_ad(i,j,k) = ltm ? 0. : st;

    });     

    Imzeli.clear(); 
    Ipzeli.clear(); 
}
