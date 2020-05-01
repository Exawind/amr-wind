#include <incflo.H>
#include <random>
#include <AMReX_Random.H>
#include "trig_ops.H"

#include "ABL.H"

using namespace amrex;

void incflo::prob_init_fluid (int lev)
{
    Box const& domain = geom[lev].Domain();
    auto const& dx = geom[lev].CellSizeArray();
    auto const& problo = geom[lev].ProbLoArray();
    auto const& probhi = geom[lev].ProbHiArray();

    auto& vel = velocity()(lev);
    auto& rho = density()(lev);
    auto& trac = temperature()(lev);
    auto& pres = pressure()(lev);

    pressure()(lev).setVal(0.0);
    grad_p()(lev).setVal(0.0);

    rho.setVal(m_ro_0);
    density().state(amr_wind::FieldState::Old)(lev).setVal(m_ro_0);

    vel.setVal(m_ic_u, 0, 1);
    vel.setVal(m_ic_v, 1, 1);
    vel.setVal(m_ic_w, 2, 1);
    trac.setVal(0.0);

    // FIXME: Ongoing refactor handle ABL/wind physics through physics interface
    if (m_probtype == 35) return;

    for (MFIter mfi(rho); mfi.isValid(); ++mfi)
    {
        const Box& vbx = mfi.validbox();
        const Box& gbx = mfi.fabbox();
        if (0 == m_probtype)
        { }
        else if (31  == m_probtype)
        {
            init_plane_poiseuille(vbx, gbx,
                                  pres.array(mfi),
                                  vel.array(mfi),
                                  rho.array(mfi),
                                  trac.array(mfi),
                                  domain, dx, problo, probhi);
        } else {
            amrex::Abort("prob_init_fluid: unknown m_probtype");
        };
    }
}

void incflo::init_plane_poiseuille (Box const& vbx, Box const& /* gbx */,
                                    Array4<Real> const& /* p */,
                                    Array4<Real> const& vel,
                                    Array4<Real> const& /* density */,
                                    Array4<Real> const& tracer,
                                    Box const& domain,
                                    GpuArray<Real, AMREX_SPACEDIM> const& /* dx */,
                                    GpuArray<Real, AMREX_SPACEDIM> const& /* problo */,
                                    GpuArray<Real, AMREX_SPACEDIM> const& /* probhi */)
{
    Real dyinv = 1.0 / domain.length(1);
    const auto dhi = amrex::ubound(domain);

    Real u = m_ic_u;
    amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        Real y = (j+0.5)*dyinv;
        vel(i,j,k,0) = 6. * u * y * (1.-y);
        vel(i,j,k,1) = 0.0;
        vel(i,j,k,2) = 0.0;

        if (i <= dhi.x/8)
            tracer(i,j,k) = 1.0;
        else
            tracer(i,j,k) = 0.0;

    });


}
