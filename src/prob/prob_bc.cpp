#include <incflo.H>

using namespace amrex;

void incflo::prob_set_inflow_velocity (int grid_id, Orientation ori, Box const& bx,
                                       Array4<Real> const& vel, int lev, Real time)
{
    if (32 == m_probtype)
    {
        Real dzinv = 1.0 / Geom(lev).Domain().length(2);
        Real v = m_ic_v;
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real z = (k+0.5)*dzinv;
            vel(i,j,k,1) = 6. * v * z * (1.-z);
        });
    }
    else
    {
        amrex::Abort("prob_set_inflow_velocity: unknown probtype");
    };
}
