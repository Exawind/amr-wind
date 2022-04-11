#include "amr-wind/incflo.H"

#include "AMReX_ParmParse.H"

using namespace amrex;

void incflo::set_background_pressure()
{
    const auto problo = geom[0].ProbLoArray();
    // determine probhi based on if mesh is stretched
    amrex::Vector<amrex::Real> probhi_physical{{0.0, 0.0, 0.0}};
    {
        amrex::ParmParse pp("geometry");
        if (pp.contains("prob_hi_physical")) {
            pp.getarr("prob_hi_physical", probhi_physical);
        } else {
            for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                probhi_physical[d] = geom[0].ProbHiArray()[d];
            }
        }
    }
    GpuArray<Real, AMREX_SPACEDIM> problen{
        {probhi_physical[0] - problo[0], probhi_physical[1] - problo[1],
         probhi_physical[2] - problo[2]}};

    amrex::Vector<amrex::Real> m_gp0{{0.0, 0.0, 0.0}};

    int delp_dir = -1;
    // (2) pressure inflow and pressure outflow
    const auto& bctype = pressure().bc_type();
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        if ((bctype[Orientation(dir, Orientation::low)] ==
                 BC::pressure_inflow and
             bctype[Orientation(dir, Orientation::high)] ==
                 BC::pressure_outflow) or
            (bctype[Orientation(dir, Orientation::high)] ==
                 BC::pressure_inflow and
             bctype[Orientation(dir, Orientation::low)] ==
                 BC::pressure_outflow)) {
            // fixme
            amrex::Abort(
                "m_gp0 is being filled with (pr_o-pr_i)/L but not used in the "
                "gradp forcing need to fix this");

            if (delp_dir == -1) {
                delp_dir = dir;
                m_gp0[dir] =
                    (pressure()
                         .bc_values()[Orientation(dir, Orientation::high)][0] -
                     pressure()
                         .bc_values()[Orientation(dir, Orientation::low)][0]) /
                    problen[dir];
            } else {
                amrex::Abort("set_background_pressure: how did this happen?");
            }
        }
    }
}
