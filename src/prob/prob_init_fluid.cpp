#include <incflo.H>

using namespace amrex;

void incflo::prob_init_fluid (int lev)
{
    auto& ld = *m_leveldata[lev];
    Box const& domain = geom[lev].Domain();
    auto const& dx = geom[lev].CellSizeArray();
    auto const& problo = geom[lev].ProbLoArray();
    auto const& probhi = geom[lev].ProbHiArray();

    ld.p.setVal(0.0);
    ld.gp.setVal(0.0);

    if (m_probtype == 0) {
        ld.density.setVal(m_ro_0);
        ld.velocity.setVal(m_ic_u, 0, 1);
        ld.velocity.setVal(m_ic_v, 1, 1);
        ld.velocity.setVal(m_ic_w, 2, 1);
        if (m_ntrac > 0) ld.tracer.setVal(0.0);
    }

    for (MFIter mfi(ld.density); mfi.isValid(); ++mfi)
    {
        const Box& vbx = mfi.validbox();
        const Box& gbx = mfi.fabbox();
        if (0 == m_probtype)
        { }
        else if (1 == m_probtype)
        {
            init_taylor_green(vbx, gbx,
                              ld.p.array(mfi),
                              ld.velocity.array(mfi),
                              ld.density.array(mfi),
                              ld.tracer.array(mfi),
                              domain, dx, problo, probhi);
        }
        else if (21 == m_probtype or 22 == m_probtype or 23 == m_probtype)
        {
            init_double_shear_layer(vbx, gbx,
                                    ld.p.array(mfi),
                                    ld.velocity.array(mfi),
                                    ld.density.array(mfi),
                                    ld.tracer.array(mfi),
                                    domain, dx, problo, probhi);
        }
        else if (32 == m_probtype)
        {
            init_plane_poiseuille(vbx, gbx,
                                  ld.p.array(mfi),
                                  ld.velocity.array(mfi),
                                  ld.density.array(mfi),
                                  ld.tracer.array(mfi),
                                  domain, dx, problo, probhi);
        }
        else
        {
            amrex::Abort("prob_init_fluid: unknown m_probtype");
        };
    }
}

void incflo::init_taylor_green (Box const& vbx, Box const& gbx,
                                Array4<Real> const& p,
                                Array4<Real> const& vel,
                                Array4<Real> const& density,
                                Array4<Real> const& tracer,
                                Box const& domain,
                                GpuArray<Real, AMREX_SPACEDIM> const& dx,
                                GpuArray<Real, AMREX_SPACEDIM> const& problo,
                                GpuArray<Real, AMREX_SPACEDIM> const& probhi)
{
    Real lrho = m_ro_0;
    amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        Real x = (i+0.5)*dx[0];
        Real y = (j+0.5)*dx[1];
        constexpr Real twopi = 2.*3.1415926535897932;
        vel(i,j,k,0) =  std::sin(twopi*x) * std::cos(twopi*y);
        vel(i,j,k,1) = -std::cos(twopi*x) * std::sin(twopi*y);
        vel(i,j,k,2) = 0.0;

        density(i,j,k) = lrho;

        const int nt = tracer.nComp();
        for (int n = 0; n < nt; ++n) {
            tracer(i,j,k,n) = 0.0;
        }
    });
}

void incflo::init_plane_poiseuille (Box const& vbx, Box const& gbx,
                                    Array4<Real> const& p,
                                    Array4<Real> const& vel,
                                    Array4<Real> const& density,
                                    Array4<Real> const& tracer,
                                    Box const& domain,
                                    GpuArray<Real, AMREX_SPACEDIM> const& dx,
                                    GpuArray<Real, AMREX_SPACEDIM> const& problo,
                                    GpuArray<Real, AMREX_SPACEDIM> const& probhi)
{
    Real dzinv = 1.0 / domain.length(2);
    const auto dlo = amrex::lbound(domain);
    const auto dhi = amrex::ubound(domain);
    Real lrho = m_ro_0;
    if (32 == m_probtype)
    {
        Real v = m_ic_v;
        amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real z = (k+0.5)*dzinv;
            vel(i,j,k,0) = 0.0;
            vel(i,j,k,1) = 6. * v * z * (1.-z);
            vel(i,j,k,2) = 0.0;

            density(i,j,k) = lrho;

            const int nt = tracer.nComp();
            for (int n = 0; n < nt; ++n) {
                tracer(i,j,k,n) = 0.0;
            }
            if (nt > 0 and j <= dhi.y/8)   tracer(i,j,k,0) = 1.0;
            if (nt > 1 and j <= dhi.y/2)   tracer(i,j,k,1) = 2.0;
            if (nt > 2 and j <= dhi.y*3/4) tracer(i,j,k,2) = 3.0;
        });
    }
    else
    {
        amrex::Abort("Unknown plane poiseuille m_probtype");
    };
}

void incflo::init_double_shear_layer (Box const& vbx, Box const& gbx,
                                      Array4<Real> const& p,
                                      Array4<Real> const& vel,
                                      Array4<Real> const& density,
                                      Array4<Real> const& tracer,
                                      Box const& domain,
                                      GpuArray<Real, AMREX_SPACEDIM> const& dx,
                                      GpuArray<Real, AMREX_SPACEDIM> const& problo,
                                      GpuArray<Real, AMREX_SPACEDIM> const& probhi)
{
    static constexpr Real twopi = 2.0 * 3.1415926535897932;
    Real lrho = m_ro_0;
    if (21 == m_probtype)
    {
        amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real x = (i+0.5) * dx[0];
            Real y = (j+0.5) * dx[1];
            vel(i,j,k,0) = std::tanh(30.0*(0.25-std::abs(y-0.5)));
            vel(i,j,k,1) = 0.05*std::sin(twopi*x);
            vel(i,j,k,2) = 0.0;

            density(i,j,k) = lrho;

            const int nt = tracer.nComp();
            for (int n = 0; n < nt; ++n) {
                tracer(i,j,k,n) = 0.0;
            }
        });
    }
    else if (22 == m_probtype)
    {
        amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real y = (j+0.5) * dx[1];
            Real z = (k+0.5) * dx[2];
            vel(i,j,k,1) = std::tanh(30.0*(0.25-std::abs(z-0.5)));
            vel(i,j,k,2) = 0.05*std::sin(twopi*y);
            vel(i,j,k,0) = 0.0;

            density(i,j,k) = lrho;

            const int nt = tracer.nComp();
            for (int n = 0; n < nt; ++n) {
                tracer(i,j,k,n) = 0.0;
            }
        });
    }
    else if (23 == m_probtype)
    {
        amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real x = (i+0.5) * dx[0];
            Real z = (k+0.5) * dx[2];
            vel(i,j,k,2) = std::tanh(30.0*(0.25-std::abs(x-0.5)));
            vel(i,j,k,0) = 0.05*std::sin(twopi*z);
            vel(i,j,k,1) = 0.0;

            density(i,j,k) = lrho;

            const int nt = tracer.nComp();
            for (int n = 0; n < nt; ++n) {
                tracer(i,j,k,n) = 0.0;
            }
        });
    }
    else
    {
        amrex::Abort("Unknown double shear layer m_probtype");
    };
}

