#include <incflo.H>

using namespace amrex;

void incflo::set_background_pressure ()
{
    m_p000 = m_ic_p;

    if (m_probtype == 11) {
        m_use_boussinesq = true;
    } else {
        const auto problo = geom[0].ProbLoArray();
        const auto probhi = geom[0].ProbHiArray();
        GpuArray<Real,AMREX_SPACEDIM> problen{probhi[0]-problo[0],
                                              probhi[1]-problo[1],
                                              probhi[2]-problo[2]};
        // There are 3 exclusive sources for background pressure gradient.
        // (1) incflo.delp in inputs
        int delp_dir = -1;
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            if (std::abs(m_delp[dir]) > std::numeric_limits<Real>::epsilon()) {
                if (delp_dir == -1) {
                    delp_dir = dir;
                    m_gp0[dir] = -m_delp[dir] / problen[dir];
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
                    m_gp0[dir] = (m_bc_pressure[Orientation(dir,Orientation::high)]
                                - m_bc_pressure[Orientation(dir,Orientation::low)]) / problen[dir];
                } else {
                    amrex::Abort("set_background_pressure: how did this happen?");
                }
            }
        }
        // (3) gravity
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            Real dpdx = m_gravity[dir] * m_ro_0;
            if (std::abs(dpdx) > std::numeric_limits<Real>::epsilon()) {
                if (delp_dir == -1) {
                    delp_dir = dir;
                    m_gp0[dir] = -dpdx;
                } else {
                    amrex::Abort("set_background_pressure: how did this happen?");
                }
            }
        }
    }
}
