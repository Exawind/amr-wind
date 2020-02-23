#include <incflo.H>
#include <AMReX_Random.H>

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

    ld.density.setVal(m_ro_0);
    ld.velocity.setVal(m_ic_u, 0, 1);
    ld.velocity.setVal(m_ic_v, 1, 1);
    ld.velocity.setVal(m_ic_w, 2, 1);
    if (m_ntrac > 0) ld.tracer.setVal(0.0);

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
        else if (3 == m_probtype)
        {
            init_taylor_green3d(vbx, gbx,
                                ld.p.array(mfi),
                                ld.velocity.array(mfi),
                                ld.density.array(mfi),
                                ld.tracer.array(mfi),
                                domain, dx, problo, probhi);
        }
        else if (4 == m_probtype)
        {
            init_couette(vbx, gbx,
                         ld.p.array(mfi),
                         ld.velocity.array(mfi),
                         ld.density.array(mfi),
                         ld.tracer.array(mfi),
                         domain, dx, problo, probhi);
        }
        else if (5 == m_probtype)
        {
            init_rayleigh_taylor(vbx, gbx,
                                 ld.p.array(mfi),
                                 ld.velocity.array(mfi),
                                 ld.density.array(mfi),
                                 ld.tracer.array(mfi),
                                 domain, dx, problo, probhi);
        }
        else if (11 == m_probtype)
        {
            init_tuscan(vbx, gbx,
                        ld.p.array(mfi),
                        ld.velocity.array(mfi),
                        ld.density.array(mfi),
                        ld.tracer.array(mfi),
                        domain, dx, problo, probhi);
        }
        else if (12 == m_probtype)
        {
            init_periodic_tracer(vbx, gbx,
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
        else if (31  == m_probtype or 32  == m_probtype or 33  == m_probtype or 
                 311 == m_probtype or 322 == m_probtype or 333 == m_probtype)
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
    amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        Real x = (i+0.5)*dx[0];
        Real y = (j+0.5)*dx[1];
        constexpr Real twopi = 2.*3.1415926535897932;
        vel(i,j,k,0) =  std::sin(twopi*x) * std::cos(twopi*y);
        vel(i,j,k,1) = -std::cos(twopi*x) * std::sin(twopi*y);
        vel(i,j,k,2) = 0.0;
    });
}

void incflo::init_taylor_green3d (Box const& vbx, Box const& gbx,
                                  Array4<Real> const& p,
                                  Array4<Real> const& vel,
                                  Array4<Real> const& density,
                                  Array4<Real> const& tracer,
                                  Box const& domain,
                                  GpuArray<Real, AMREX_SPACEDIM> const& dx,
                                  GpuArray<Real, AMREX_SPACEDIM> const& problo,
                                  GpuArray<Real, AMREX_SPACEDIM> const& probhi)
{
    amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        Real x = (i+0.5)*dx[0];
        Real y = (j+0.5)*dx[1];
        Real z = (k+0.5)*dx[2];
        constexpr Real twopi = 2.*3.1415926535897932;
        vel(i,j,k,0) =  std::sin(twopi*x) * std::cos(twopi*y) * cos(twopi*z);
        vel(i,j,k,1) = -std::cos(twopi*x) * std::sin(twopi*y) * cos(twopi*z);
        vel(i,j,k,2) = 0.0;
    });
}

void incflo::init_couette (Box const& vbx, Box const& gbx,
                           Array4<Real> const& p,
                           Array4<Real> const& vel,
                           Array4<Real> const& density,
                           Array4<Real> const& tracer,
                           Box const& domain,
                           GpuArray<Real, AMREX_SPACEDIM> const& dx,
                           GpuArray<Real, AMREX_SPACEDIM> const& problo,
                           GpuArray<Real, AMREX_SPACEDIM> const& probhi)
{
    Real num_cells_y = static_cast<Real>(domain.length(1));
    amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        Real y = (j+0.5) / num_cells_y;
        vel(i,j,k,0) *= (y-0.5);
        vel(i,j,k,1) = 0.0;
        vel(i,j,k,2) = 0.0;
    });
}

void incflo::init_rayleigh_taylor (Box const& vbx, Box const& gbx,
                                   Array4<Real> const& p,
                                   Array4<Real> const& vel,
                                   Array4<Real> const& density,
                                   Array4<Real> const& tracer,
                                   Box const& domain,
                                   GpuArray<Real, AMREX_SPACEDIM> const& dx,
                                   GpuArray<Real, AMREX_SPACEDIM> const& problo,
                                   GpuArray<Real, AMREX_SPACEDIM> const& probhi)
{
    static constexpr Real pi = 3.1415926535897932;
    static constexpr Real rho_1 = 0.5;
    static constexpr Real rho_2 = 2.0;
    const Real splitx = 0.5*(problo[0] + probhi[0]);
    const Real splity = 0.5*(problo[1] + probhi[1]);
    const Real L_x = probhi[0] - problo[0];

    amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        Real x = problo[0] + (i+0.5)*dx[0];
        Real y = problo[1] + (j+0.5)*dx[1];
        Real z = problo[2] + (k+0.5)*dx[2];
        const Real r2d = amrex::min(std::hypot((x-splitx),(y-splity)), 0.5*L_x);
        const Real pertheight = 0.5 - 0.01*std::cos(2.0*pi*r2d/L_x);
        density(i,j,k) = rho_1 + ((rho_2-rho_1)/2.0)*(1.0+std::tanh((z-pertheight)/0.005));
        vel(i,j,k,0) = 0.0;
        vel(i,j,k,1) = 0.0;
        vel(i,j,k,2) = 0.0;
    });
}

void incflo::init_tuscan (Box const& vbx, Box const& gbx,
                          Array4<Real> const& p,
                          Array4<Real> const& vel,
                          Array4<Real> const& density,
                          Array4<Real> const& tracer,
                          Box const& domain,
                          GpuArray<Real, AMREX_SPACEDIM> const& dx,
                          GpuArray<Real, AMREX_SPACEDIM> const& problo,
                          GpuArray<Real, AMREX_SPACEDIM> const& probhi)
{
    int half_num_cells = domain.length(2) / 2;
    amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        vel(i,j,k,0) = 0.0;
        vel(i,j,k,1) = 0.0;
        vel(i,j,k,2) = 0.0;
        density(i,j,k) = 1.0;
        if (k <= half_num_cells) {
            tracer(i,j,k) = 0.0;
        } else {
            tracer(i,j,k) = 0.01;
        }
    });
}
void incflo::init_periodic_tracer (Box const& vbx, Box const& gbx,
                                   Array4<Real> const& p,
                                   Array4<Real> const& vel,
                                   Array4<Real> const& density,
                                   Array4<Real> const& tracer,
                                   Box const& domain,
                                   GpuArray<Real, AMREX_SPACEDIM> const& dx,
                                   GpuArray<Real, AMREX_SPACEDIM> const& problo,
                                   GpuArray<Real, AMREX_SPACEDIM> const& probhi)
{
    Real L = probhi[0]-problo[0];
    Real C = 2.*3.1415926535897932 / L;
    amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        constexpr Real A = 1.0;
        Real x = (i+0.5)*dx[0];
        Real y = (j+0.5)*dx[1];
        Real z = (k+0.5)*dx[2];
        vel(i,j,k,0) = 1.0;
        vel(i,j,k,1) = 0.1*(std::sin(C*(x+z) - 0.00042) + 1.0) * std::exp(y);
        vel(i,j,k,2) = 0.1*(std::sin(C*(x+y) - 0.00042) + 1.0) * std::exp(z);
        tracer(i,j,k) = A *(std::sin(C*(y+z) - 0.00042) + 1.0) * std::exp(x);
    });
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
     if (21 == m_probtype)
    {
        amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real x = (i+0.5) * dx[0];
            Real y = (j+0.5) * dx[1];
            vel(i,j,k,0) = std::tanh(30.0*(0.25-std::abs(y-0.5)));
            vel(i,j,k,1) = 0.05*std::sin(twopi*x);
            vel(i,j,k,2) = 0.0;
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
        });
    }
    else
    {
        amrex::Abort("Unknown double shear layer m_probtype");
    };
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
    Real dxinv = 1.0 / domain.length(0);
    Real dyinv = 1.0 / domain.length(1);
    Real dzinv = 1.0 / domain.length(2);
    const auto dlo = amrex::lbound(domain);
    const auto dhi = amrex::ubound(domain);

    if (31 == m_probtype)
    {
        Real u = m_ic_u;
        amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real y = (j+0.5)*dyinv;
            vel(i,j,k,0) = 6. * u * y * (1.-y);
            vel(i,j,k,1) = 0.0;
            vel(i,j,k,2) = 0.0;

            const int nt = tracer.nComp();
            for (int n = 0; n < nt; ++n) {
                tracer(i,j,k,n) = 0.0;
            }
            if (nt > 0 and i <= dhi.x/8)   tracer(i,j,k,0) = 1.0;
            if (nt > 1 and i <= dhi.x/2)   tracer(i,j,k,1) = 2.0;
            if (nt > 2 and i <= dhi.x*3/4) tracer(i,j,k,2) = 3.0;
        });
    }
    else if (311 == m_probtype)
    {
        Real u = m_ic_u;
        amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real z = (k+0.5)*dzinv;
            vel(i,j,k,0) = 6. * u * z * (1.-z);
            vel(i,j,k,1) = 0.0;
            vel(i,j,k,2) = 0.0;

            const int nt = tracer.nComp();
            for (int n = 0; n < nt; ++n) {
                tracer(i,j,k,n) = 0.0;
            }
            if (nt > 0 and i <= dhi.x/8)   tracer(i,j,k,0) = 1.0;
            if (nt > 1 and i <= dhi.x/2)   tracer(i,j,k,1) = 2.0;
            if (nt > 2 and i <= dhi.x*3/4) tracer(i,j,k,2) = 3.0;
        });
    }
    else if (32 == m_probtype)
    {
        Real v = m_ic_v;
        amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real z = (k+0.5)*dzinv;
            vel(i,j,k,0) = 0.0;
            vel(i,j,k,1) = 6. * v * z * (1.-z);
            vel(i,j,k,2) = 0.0;

            const int nt = tracer.nComp();
            for (int n = 0; n < nt; ++n) {
                tracer(i,j,k,n) = 0.0;
            }
            if (nt > 0 and j <= dhi.y/8)   tracer(i,j,k,0) = 1.0;
            if (nt > 1 and j <= dhi.y/2)   tracer(i,j,k,1) = 2.0;
            if (nt > 2 and j <= dhi.y*3/4) tracer(i,j,k,2) = 3.0;
        });
    }
    else if (322 == m_probtype)
    {
        Real v = m_ic_v;
        amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real x = (i+0.5)*dxinv;
            vel(i,j,k,0) = 0.0;
            vel(i,j,k,1) = 6. * v * x * (1.-x);
            vel(i,j,k,2) = 0.0;

            const int nt = tracer.nComp();
            for (int n = 0; n < nt; ++n) {
                tracer(i,j,k,n) = 0.0;
            }
            if (nt > 0 and j <= dhi.y/8)   tracer(i,j,k,0) = 1.0;
            if (nt > 1 and j <= dhi.y/2)   tracer(i,j,k,1) = 2.0;
            if (nt > 2 and j <= dhi.y*3/4) tracer(i,j,k,2) = 3.0;
        });
    }
    else if (33 == m_probtype)
    {
        Real w = m_ic_w;
        amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real x = (i+0.5)*dxinv;
            vel(i,j,k,0) = 0.0;
            vel(i,j,k,1) = 0.0;
            vel(i,j,k,2) = 6. * w * x * (1.-x);

            const int nt = tracer.nComp();
            for (int n = 0; n < nt; ++n) {
                tracer(i,j,k,n) = 0.0;
            }
            if (nt > 0 and k <= dhi.z/8)   tracer(i,j,k,0) = 1.0;
            if (nt > 1 and k <= dhi.z/2)   tracer(i,j,k,1) = 2.0;
            if (nt > 2 and k <= dhi.z*3/4) tracer(i,j,k,2) = 3.0;
        });
    }
    else if (333 == m_probtype)
    {
        Real w = m_ic_w;
        amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real y = (j+0.5)*dyinv;
            vel(i,j,k,0) = 0.0;
            vel(i,j,k,1) = 0.0;
            vel(i,j,k,2) = 6. * w * y * (1.-y);

            const int nt = tracer.nComp();
            for (int n = 0; n < nt; ++n) {
                tracer(i,j,k,n) = 0.0;
            }
            if (nt > 0 and k <= dhi.z/8)   tracer(i,j,k,0) = 1.0;
            if (nt > 1 and k <= dhi.z/2)   tracer(i,j,k,1) = 2.0;
            if (nt > 2 and k <= dhi.z*3/4) tracer(i,j,k,2) = 3.0;
        });
    }
    else
    {
        amrex::Abort("Unknown plane poiseuille m_probtype");
    };
}
