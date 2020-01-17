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
#pragma omp parallel
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

//
//  These subroutines set the BCs for the vel_arr components only.
//

void
incflo::incflo_set_velocity_bcs (Real time,
                                 Vector< std::unique_ptr<MultiFab> > & vel_in)
{
  BL_PROFILE("incflo::incflo_set_velocity_bcs()");

  for (int lev = 0; lev <= finest_level; lev++)
  {
     // Set all values outside the domain to covered_val just to avoid use of undefined
     vel_in[lev]->setDomainBndry(covered_val,geom[lev]);

     vel_in[lev] -> FillBoundary (geom[lev].periodicity());
     Box domain(geom[lev].Domain());

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
     for (MFIter mfi(*vel_in[lev]); mfi.isValid(); ++mfi) {
         set_velocity_bcs(time, lev, (*vel_in[lev])[mfi], domain);
     }

#ifdef AMREX_USE_EB
     EB_set_covered(*vel_in[lev], 0, vel_in[lev]->nComp(), vel_in[lev]->nGrow(), covered_val);
#endif

     // Do this after as well as before to pick up terms that got updated in the call above
     vel_in[lev] -> FillBoundary (geom[lev].periodicity());
  }
}

void
incflo::set_velocity_bcs(Real time,
                         const int lev,
                         FArrayBox& vel_fab,
                         const Box& domain) const
{
}
