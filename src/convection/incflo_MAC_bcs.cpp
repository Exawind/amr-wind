#include <incflo_MAC_bcs.H>

using namespace amrex;

void incflo_set_mac_bcs (Box const& domain,
                         Box const& ubx,
                         Box const& vbx,
                         Box const& wbx,
                         Array4<Real> const& u,
                         Array4<Real> const& v,
                         Array4<Real> const& w,
                         Array4<Real const> const& vel,
                         Vector<BCRec> const& h_bcrec,
                         BCRec const* d_bcrec)
{
    int idim = 0;
    if (h_bcrec[idim].lo(idim) == BCType::ext_dir and
        domain.smallEnd(idim) == ubx.smallEnd(idim))
    {
        amrex::ParallelFor(amrex::bdryLo(ubx,idim),
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            u(i,j,k) = vel(i-1,j,k,0);
        });
    }

    if (h_bcrec[idim].hi(idim) == BCType::ext_dir and
        domain.bigEnd(idim)+1 == ubx.bigEnd(idim))
    {
        amrex::ParallelFor(amrex::bdryHi(ubx,idim),
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            u(i,j,k) = vel(i,j,k,0);
        });
    }

    idim = 1;
    if (h_bcrec[idim].lo(idim) == BCType::ext_dir and
        domain.smallEnd(idim) == vbx.smallEnd(idim))
    {
        amrex::ParallelFor(amrex::bdryLo(vbx,idim),
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            v(i,j,k) = vel(i,j-1,k,1);
        });
    }

    if (h_bcrec[idim].hi(idim) == BCType::ext_dir and
        domain.bigEnd(idim)+1 == vbx.bigEnd(idim))
    {
        amrex::ParallelFor(amrex::bdryHi(vbx,idim),
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            v(i,j,k) = vel(i,j,k,1);
        });
    }

    idim = 2;
    if (h_bcrec[idim].lo(idim) == BCType::ext_dir and
        domain.smallEnd(idim) == wbx.smallEnd(idim))
    {
        amrex::ParallelFor(amrex::bdryLo(wbx,idim),
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            w(i,j,k) = vel(i,j,k-1,2);
        });
    }

    if (h_bcrec[idim].hi(idim) == BCType::ext_dir and
        domain.bigEnd(idim)+1 == wbx.bigEnd(idim))
    {
        amrex::ParallelFor(amrex::bdryHi(wbx,idim),
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            w(i,j,k) = vel(i,j,k,2);
        });
    }
}


