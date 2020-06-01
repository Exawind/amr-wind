#include "amr-wind/incflo.H"

using namespace amrex;

void incflo::prob_set_inflow_velocity (int /* grid_id */, Orientation ori, Box const& bx,
                                       Array4<Real> const& vel, int lev, Real /* time */,
                                       Real bcv_in)
{
    if (31 == m_probtype)
    {
        Real dyinv = 1.0 / Geom(lev).Domain().length(1);
        Real u = 1.0;
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real y = (j+0.5)*dyinv;
            vel(i,j,k,0) = 6. * u * y * (1.-y);
        });
    }
    else
    {
        const int  dir = ori.coordDir();
        const Real bcv = bcv_in;
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept 
        { 
            vel(i,j,k,dir) = bcv; 
        });
    };
}
