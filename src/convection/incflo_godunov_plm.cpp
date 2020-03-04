#include <incflo_convection_K.H>
#include <incflo.H>

using namespace amrex;

void incflo::predict_plm (int lev, Box const& ubx, Box const& vbx, Box const& wbx,
                          Array4<Real> const& Imx, Array4<Real> const& Ipx,
                          Array4<Real> const& Imy, Array4<Real> const& Ipy,
                          Array4<Real> const& Imz, Array4<Real> const& Ipz,
                          Array4<Real const> const& vcc)
{
    const Real dx = geom[lev].CellSize(0);
    const Real dy = geom[lev].CellSize(1);
    const Real dz = geom[lev].CellSize(2);
    Real l_dt = m_dt;
    Real dtdx = l_dt/dx;
    Real dtdy = l_dt/dy;
    Real dtdz = l_dt/dz;

    const Box& domain_box = geom[lev].Domain();
    const int domain_ilo = domain_box.smallEnd(0);
    const int domain_ihi = domain_box.bigEnd(0);
    const int domain_jlo = domain_box.smallEnd(1);
    const int domain_jhi = domain_box.bigEnd(1);
    const int domain_klo = domain_box.smallEnd(2);
    const int domain_khi = domain_box.bigEnd(2);

    auto const bc_ilo = m_bc_type[Orientation(Direction::x,Orientation::low)];
    auto const bc_ihi = m_bc_type[Orientation(Direction::x,Orientation::high)];
    auto const bc_jlo = m_bc_type[Orientation(Direction::y,Orientation::low)];
    auto const bc_jhi = m_bc_type[Orientation(Direction::y,Orientation::high)];
    auto const bc_klo = m_bc_type[Orientation(Direction::z,Orientation::low)];
    auto const bc_khi = m_bc_type[Orientation(Direction::z,Orientation::high)];

    bool extdir_ilo = (bc_ilo == BC::mass_inflow) or (bc_ilo == BC::no_slip_wall);
    bool extdir_ihi = (bc_ihi == BC::mass_inflow) or (bc_ihi == BC::no_slip_wall);
    bool extdir_jlo = (bc_jlo == BC::mass_inflow) or (bc_jlo == BC::no_slip_wall);
    bool extdir_jhi = (bc_jhi == BC::mass_inflow) or (bc_jhi == BC::no_slip_wall);
    bool extdir_klo = (bc_klo == BC::mass_inflow) or (bc_klo == BC::no_slip_wall);
    bool extdir_khi = (bc_khi == BC::mass_inflow) or (bc_khi == BC::no_slip_wall);

    // At an ext_dir boundary, the boundary value is on the face, not cell center.

    if ((extdir_ilo and domain_ilo >= ubx.smallEnd(0)-1) or
        (extdir_ihi and domain_ihi <= ubx.bigEnd(0)))
    {
        amrex::ParallelFor(ubx, [vcc,extdir_ilo,extdir_ihi,domain_ilo,domain_ihi,Imx,Ipx,dtdx]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real upls = vcc(i  ,j,k,0) + 0.5 * (-1.0 - vcc(i  ,j,k,0) * dtdx) * 
                incflo_xslope_extdir(i,j,k,0,vcc, extdir_ilo, extdir_ihi, domain_ilo, domain_ihi);
            Real umns = vcc(i-1,j,k,0) + 0.5 * ( 1.0 - vcc(i-1,j,k,0) * dtdx) * 
                incflo_xslope_extdir(i-1,j,k,0,vcc, extdir_ilo, extdir_ihi, domain_ilo, domain_ihi);

            if (extdir_ilo and i == domain_ilo) {
                umns = vcc(i-1,j,k,0);
                upls = vcc(i-1,j,k,0);
            } else if (extdir_ihi and i == domain_ihi+1) {
                umns = vcc(i,j,k,0);
                upls = vcc(i,j,k,0);
            }

            Ipx(i-1,j,k,0) = umns;
            Imx(i  ,j,k,0) = upls;
        });
    }
    else
    {
        amrex::ParallelFor(ubx, [vcc,Ipx,Imx,dtdx]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real upls = vcc(i  ,j,k,0) + 0.5 * (-1.0 - vcc(i  ,j,k,0) * dtdx) * 
                incflo_xslope(i  ,j,k,0,vcc);
            Real umns = vcc(i-1,j,k,0) + 0.5 * ( 1.0 - vcc(i-1,j,k,0) * dtdx) * 
                incflo_xslope(i-1,j,k,0,vcc);

            Ipx(i-1,j,k,0) = umns;
            Imx(i  ,j,k,0) = upls;
        });
    }

    if ((extdir_jlo and domain_jlo >= vbx.smallEnd(1)-1) or
        (extdir_jhi and domain_jhi <= vbx.bigEnd(1)))
    {
        amrex::ParallelFor(vbx, [vcc,extdir_jlo,extdir_jhi,domain_jlo,domain_jhi,Imy,Ipy,dtdy]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real vpls = vcc(i,j  ,k,1) + 0.5 * (-1.0 - vcc(i,j  ,k,1) * dtdy) * 
                incflo_yslope_extdir(i,j,k,1,vcc, extdir_jlo, extdir_jhi, domain_jlo, domain_jhi);
            Real vmns = vcc(i,j-1,k,1) + 0.5 * ( 1.0 - vcc(i,j-1,k,1) * dtdy) * 
                incflo_yslope_extdir(i,j-1,k,1,vcc, extdir_jlo, extdir_jhi, domain_jlo, domain_jhi);

            if (extdir_jlo and j == domain_jlo) {
                vmns = vcc(i,j-1,k,1);
                vpls = vcc(i,j-1,k,1);
            } else if (extdir_jhi and j == domain_jhi+1) {
                vmns = vcc(i,j,k,1);
                vpls = vcc(i,j,k,1);
            }

            Ipy(i,j-1,k,1) = vmns;
            Imy(i,j  ,k,1) = vpls;
        });
    }
    else
    {
        amrex::ParallelFor(vbx, [vcc,Ipy,Imy,dtdy]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real vpls = vcc(i,j  ,k,1) + 0.5 * (-1.0 - vcc(i,j  ,k,1) * dtdy) * 
                incflo_yslope(i,j  ,k,1,vcc);
            Real vmns = vcc(i,j-1,k,1) + 0.5 * ( 1.0 - vcc(i,j-1,k,1) * dtdy) * 
                incflo_yslope(i,j-1,k,1,vcc);

            Ipy(i,j-1,k,1) = vmns;
            Imy(i,j  ,k,1) = vpls;
        });
    }

    if ((extdir_klo and domain_klo >= wbx.smallEnd(2)-1) or
        (extdir_khi and domain_khi <= wbx.bigEnd(2)))
    {
        amrex::ParallelFor(wbx, [vcc,extdir_klo,extdir_khi,domain_klo,domain_khi,Ipz,Imz,dtdz]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real wpls = vcc(i,j,k  ,2) + 0.5 * (-1.0 - vcc(i,j,k  ,2) * dtdz) * 
                incflo_zslope_extdir(i,j,k,2,vcc, extdir_klo, extdir_khi, domain_klo, domain_khi);
            Real wmns = vcc(i,j,k-1,2) + 0.5 * ( 1.0 - vcc(i,j,k-1,2) * dtdz) * 
                incflo_zslope_extdir(i,j,k-1,2,vcc, extdir_klo, extdir_khi, domain_klo, domain_khi);

            if (extdir_klo and k == domain_klo) {
                wmns = vcc(i,j,k-1,2);
                wpls = vcc(i,j,k-1,2);
            } else if (extdir_khi and k == domain_khi+1) {
                wmns = vcc(i,j,k,2);
                wpls = vcc(i,j,k,2);
            }

            Ipz(i,j,k-1,2) = wmns;
            Imz(i,j,k  ,2) = wpls;
        });
    }
    else
    {
        amrex::ParallelFor(wbx, [vcc,Ipz,Imz,dtdz]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real wpls = vcc(i,j,k  ,2) + 0.5 * (-1.0 - vcc(i,j,k  ,2) * dtdz) * 
                incflo_zslope(i,j,k  ,2,vcc);
            Real wmns = vcc(i,j,k-1,2) + 0.5 * ( 1.0 - vcc(i,j,k-1,2) * dtdz) * 
                incflo_zslope(i,j,k-1,2,vcc);

            Ipz(i,j,k-1,2) = wmns;
            Imz(i,j,k  ,2) = wpls;
        });
    }
}

