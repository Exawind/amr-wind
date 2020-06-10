#include "amr-wind/convection/incflo_convection_K.H"
#include "amr-wind/convection/MOL.H"
#include "amr-wind/utilities/bc_ops.H"

using namespace amrex;


void mol::predict_vels_on_faces (int lev, Box const& ubx, Box const& vbx, Box const& wbx,
                                 Array4<Real> const& u, Array4<Real> const& v, Array4<Real> const& w,
                                 Array4<Real const> const& vcc,
                                 Vector<BCRec> const& h_bcrec,
                                 BCRec const* d_bcrec,
                                 Vector<Geometry> geom)
{
    BL_PROFILE("amr-wind::mol::predict_vels_on_faces");

    constexpr Real small_vel = 1.e-10;

    const int ncomp = AMREX_SPACEDIM; // This is only used because h_bcrec and d_bcrec hold the
                                      // bc's for all three velocity components

    const Box& domain_box = geom[lev].Domain();
    const int domain_ilo = domain_box.smallEnd(0);
    const int domain_ihi = domain_box.bigEnd(0);
    const int domain_jlo = domain_box.smallEnd(1);
    const int domain_jhi = domain_box.bigEnd(1);
    const int domain_klo = domain_box.smallEnd(2);
    const int domain_khi = domain_box.bigEnd(2);

    // At an ext_dir boundary, the boundary value is on the face, not cell center.

    auto extdir_lohi = amr_wind::utils::has_extdir(h_bcrec.data(), ncomp, static_cast<int>(Direction::x));
    bool has_extdir_lo = extdir_lohi.first;
    bool has_extdir_hi = extdir_lohi.second;

    if ((has_extdir_lo and domain_ilo >= ubx.smallEnd(0)-1) or
        (has_extdir_hi and domain_ihi <= ubx.bigEnd(0)))
    {
        amrex::ParallelFor(ubx, [vcc,domain_ilo,domain_ihi,u,d_bcrec]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            bool extdir_ilo = d_bcrec[0].lo(0) == BCType::ext_dir;
            bool extdir_ihi = d_bcrec[0].hi(0) == BCType::ext_dir;

            Real upls = vcc(i,j,k,0) - 0.5 * incflo_xslope_extdir
                (i,j,k,0,vcc, extdir_ilo, extdir_ihi, domain_ilo, domain_ihi);
            Real umns = vcc(i-1,j,k,0) + 0.5 * incflo_xslope_extdir
                (i-1,j,k,0,vcc, extdir_ilo, extdir_ihi, domain_ilo, domain_ihi);
            if (umns < 0.0 and upls > 0.0) {
                u(i,j,k) = 0.0;
            } else {
                Real avg = 0.5 * (upls + umns);
                if (amrex::Math::abs(avg) < small_vel) {
                    u(i,j,k) = 0.0;
                } else if (avg > 0.0) {
                    u(i,j,k) = umns;
                } else {
                    u(i,j,k) = upls;
                }
            }

            if (extdir_ilo and i == domain_ilo) {
                u(i,j,k) = vcc(i-1,j,k,0);
            } else if (extdir_ihi and i == domain_ihi+1) {
                u(i,j,k) = vcc(i,j,k,0);
            }
        });
    }
    else
    {
        amrex::ParallelFor(ubx, [vcc,u]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real upls = vcc(i  ,j,k,0) - 0.5 * incflo_xslope(i  ,j,k,0,vcc);
            Real umns = vcc(i-1,j,k,0) + 0.5 * incflo_xslope(i-1,j,k,0,vcc);
            if (umns < 0.0 and upls > 0.0) {
                u(i,j,k) = 0.0;
            } else {
                Real avg = 0.5 * (upls + umns);
                if (amrex::Math::abs(avg) < small_vel) {
                    u(i,j,k) = 0.0;
                } else if (avg > 0.0) {
                    u(i,j,k) = umns;
                } else {
                    u(i,j,k) = upls;
                }
            }
        });
    }

    extdir_lohi = amr_wind::utils::has_extdir(h_bcrec.data(), ncomp, static_cast<int>(Direction::y));
    has_extdir_lo = extdir_lohi.first;
    has_extdir_hi = extdir_lohi.second;

    if ((has_extdir_lo and domain_jlo >= vbx.smallEnd(1)-1) or
        (has_extdir_hi and domain_jhi <= vbx.bigEnd(1)))
    {
        amrex::ParallelFor(vbx, [vcc,domain_jlo,domain_jhi,v,d_bcrec]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            bool extdir_jlo = d_bcrec[1].lo(1) == BCType::ext_dir;
            bool extdir_jhi = d_bcrec[1].hi(1) == BCType::ext_dir;

            Real vpls = vcc(i,j,k,1) - 0.5 * incflo_yslope_extdir
                (i,j,k,1,vcc, extdir_jlo, extdir_jhi, domain_jlo, domain_jhi);
            Real vmns = vcc(i,j-1,k,1) + 0.5 * incflo_yslope_extdir
                (i,j-1,k,1,vcc, extdir_jlo, extdir_jhi, domain_jlo, domain_jhi);
            if (vmns < 0.0 and vpls > 0.0) {
                v(i,j,k) = 0.0;
            } else {
                Real avg = 0.5 * (vpls + vmns);
                if (amrex::Math::abs(avg) < small_vel) {
                    v(i,j,k) = 0.0;
                } else if (avg > 0.0) {
                    v(i,j,k) = vmns;
                } else {
                    v(i,j,k) = vpls;
                }
            }

            if (extdir_jlo and j == domain_jlo) {
                v(i,j,k) = vcc(i,j-1,k,1);
            } else if (extdir_jhi and j == domain_jhi+1) {
                v(i,j,k) = vcc(i,j,k,1);
            }
        });
    }
    else
    {
        amrex::ParallelFor(vbx, [vcc,v]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real vpls = vcc(i,j  ,k,1) - 0.5 * incflo_yslope(i,j  ,k,1,vcc);
            Real vmns = vcc(i,j-1,k,1) + 0.5 * incflo_yslope(i,j-1,k,1,vcc);
            if (vmns < 0.0 and vpls > 0.0) {
                v(i,j,k) = 0.0;
            } else {
                Real avg = 0.5 * (vpls + vmns);
                if (amrex::Math::abs(avg) < small_vel) {
                    v(i,j,k) = 0.0;
                } else if (avg > 0.0) {
                    v(i,j,k) = vmns;
                } else {
                    v(i,j,k) = vpls;
                }
            }
        });
    }

    extdir_lohi = amr_wind::utils::has_extdir(h_bcrec.data(), ncomp, static_cast<int>(Direction::z));
    has_extdir_lo = extdir_lohi.first;
    has_extdir_hi = extdir_lohi.second;

    if ((has_extdir_lo and domain_klo >= wbx.smallEnd(2)-1) or
        (has_extdir_hi and domain_khi <= wbx.bigEnd(2)))
    {
        amrex::ParallelFor(wbx, [vcc,domain_klo,domain_khi,w,d_bcrec]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            bool extdir_klo = d_bcrec[2].lo(2) == BCType::ext_dir;
            bool extdir_khi = d_bcrec[2].hi(2) == BCType::ext_dir;

            Real wpls = vcc(i,j,k,2) - 0.5 * incflo_zslope_extdir
                (i,j,k,2,vcc, extdir_klo, extdir_khi, domain_klo, domain_khi);
            Real wmns = vcc(i,j,k-1,2) + 0.5 * incflo_zslope_extdir(
                i,j,k-1,2,vcc, extdir_klo, extdir_khi, domain_klo, domain_khi);
            if (wmns < 0.0 and wpls > 0.0) {
                w(i,j,k) = 0.0;
            } else {
                Real avg = 0.5 * (wpls + wmns);
                if (amrex::Math::abs(avg) < small_vel) {
                    w(i,j,k) = 0.0;
                } else if (avg > 0.0) {
                    w(i,j,k) = wmns;
                } else {
                    w(i,j,k) = wpls;
                }
            }

            if (extdir_klo and k == domain_klo) {
                w(i,j,k) = vcc(i,j,k-1,2);
            } else if (extdir_khi and k == domain_khi+1) {
                w(i,j,k) = vcc(i,j,k,2);
            }
        });
    }
    else
    {
        amrex::ParallelFor(wbx, [vcc,w]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real wpls = vcc(i,j,k  ,2) - 0.5 * incflo_zslope(i,j,k  ,2,vcc);
            Real wmns = vcc(i,j,k-1,2) + 0.5 * incflo_zslope(i,j,k-1,2,vcc);
            if (wmns < 0.0 and wpls > 0.0) {
                w(i,j,k) = 0.0;
            } else {
                Real avg = 0.5 * (wpls + wmns);
                if (amrex::Math::abs(avg) < small_vel) {
                    w(i,j,k) = 0.0;
            } else if (avg > 0.0) {
                    w(i,j,k) = wmns;
                } else {
                    w(i,j,k) = wpls;
                }
            }
        });
    }
}
