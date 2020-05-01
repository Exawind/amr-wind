#include <incflo.H>

using namespace amrex;

void
incflo::set_inflow_velocity (int lev, amrex::Real time, MultiFab& vel, int nghost)
{
    auto& velocity = icns().fields().field;
    auto& bctype = velocity.bc_type();
    Geometry const& gm = Geom(lev);
    Box const& domain = gm.growPeriodicDomain(nghost);
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        Orientation olo(dir,Orientation::low);
        Orientation ohi(dir,Orientation::high);
        if (bctype[olo] == BC::mass_inflow or bctype[ohi] == BC::mass_inflow) {
            Box dlo = (bctype[olo] == BC::mass_inflow) ? amrex::adjCellLo(domain,dir,nghost) : Box();
            Box dhi = (bctype[ohi] == BC::mass_inflow) ? amrex::adjCellHi(domain,dir,nghost) : Box();
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
                    const Real bcv = velocity.bc_values()[olo][dir];
                    prob_set_inflow_velocity(gid, olo, blo, v, lev, time, bcv);
                }
                if (bhi.ok()) {
                    const Real bcv = velocity.bc_values()[ohi][dir];
                    prob_set_inflow_velocity(gid, ohi, bhi, v, lev, time, bcv);
                }
            }
        }
    }
}
