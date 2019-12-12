#include <incflo.H>

using namespace amrex;

void incflo::set_background_pressure (int lev)
{
    auto& ld = *m_leveldata[lev];

    p000 = ic_p;

    if (probtype == 11) {
        use_boussinesq = true;
    } else {
        const auto problo = geom[lev].ProbLoArray();
        const auto probhi = geom[lev].ProbHiArray();
        GpuArray<Real,AMREX_SPACEDIM> problen{probhi[0]-problo[0],
                                              probhi[1]-problo[1],
                                              probhi[2]-problo[2]};
        // There are 3 exclusive sources for background pressure gradient.
        // (1) incflo.delp in inputs
        int delp_dir = -1;
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            if (std::abs(delp[dir]) > std::numeric_limits<Real>::epsilon()) {
                if (delp_dir == -1) {
                    delp_dir = dir;
                    gp0[dir] = delp[dir] / problen[dir];
                } else {
                    amrex::Abort("set_background_pressure: how did this happen?");
                }
            }
        }
        // (2) pressure inflow and pressure outflow
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            if ((m_bc_type[Orientation(dir,Orientation::low)] == BC::pressure_inflow and
                 m_bc_type[Orientation(dir,Orientation::high)] == BC::pressure_outflow) or
                (m_bc_type[Orientation(dir,Orientation::high)] == BC::pressure_inflow and
                 m_bc_type[Orientation(dir,Orientation::low)] == BC::pressure_outflow))
            {
                if (delp_dir == -1) {
                    delp_dir = dir;
                    gp0[dir] = (m_bc_pressure[Orientation(dir,Orientation::high)]
                                - m_bc_pressure[Orientation(dir,Orientation::low)]) / problen[dir];
                } else {
                    amrex::Abort("set_background_pressure: how did this happen?");
                }
            }
        }
        // (3) gravity
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            Real dpdx = gravity[dir] * ro_0;
            if (std::abs(dpdx) > std::numeric_limits<Real>::epsilon()) {
                if (delp_dir == -1) {
                    delp_dir = dir;
                    gp0[dir] = dpdx;
                } else {
                    amrex::Abort("set_background_pressure: how did this happen?");
                }
            }
        }
    }
}
