#ifdef AMREX_USE_EB
#include <AMReX_EBMultiFabUtil.H>
#endif

#include <incflo.H>

using namespace amrex;

void
incflo::set_inflow_velocity (int lev, amrex::Real time, MultiFab& vel, int nghost)
{
    Geometry const& gm = Geom(lev);
    Box const& domain = gm.growPeriodicDomain(nghost);
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        Orientation olo(dir,Orientation::low);
        Orientation ohi(dir,Orientation::high);
        if (m_bc_type[olo] == BC::mass_inflow or m_bc_type[ohi] == BC::mass_inflow) {
            Box dlo = (m_bc_type[olo] == BC::mass_inflow) ? amrex::adjCellLo(domain,dir,nghost) : Box();
            Box dhi = (m_bc_type[ohi] == BC::mass_inflow) ? amrex::adjCellHi(domain,dir,nghost) : Box();
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(vel); mfi.isValid(); ++mfi) {
                Box const& gbx = amrex::grow(mfi.validbox(),nghost);
                Box blo = gbx & dlo;
                Box bhi = gbx & dhi;
                Array4<Real> const& v = vel[mfi].array();
                int gid = mfi.index();
                if (blo.ok()) {
                    prob_set_inflow_velocity(gid, olo, blo, v, lev, time);
                }
                if (bhi.ok()) {
                    prob_set_inflow_velocity(gid, ohi, bhi, v, lev, time);
                }
            }
        }
    }
}
