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
    auto& trac = tracer()(lev);
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
    if (m_probtype == 35 || m_probtype == 11) return;

    for (MFIter mfi(rho); mfi.isValid(); ++mfi)
    {
        const Box& vbx = mfi.validbox();
        const Box& gbx = mfi.fabbox();
        if (0 == m_probtype)
        { }
        else if (1 == m_probtype)
        {
            init_taylor_green(vbx, gbx,
                              pres.array(mfi),
                              vel.array(mfi),
                              rho.array(mfi),
                              trac.array(mfi),
                              domain, dx, problo, probhi);
        }
        else if (3 == m_probtype)
        {
            init_taylor_green3d(vbx, gbx,
                                pres.array(mfi),
                                vel.array(mfi),
                                rho.array(mfi),
                                trac.array(mfi),
                                domain, dx, problo, probhi);
        }
        else if (4 == m_probtype)
        {
            init_couette(vbx, gbx,
                         pres.array(mfi),
                         vel.array(mfi),
                         rho.array(mfi),
                         trac.array(mfi),
                         domain, dx, problo, probhi);
        }
        else if (5 == m_probtype)
        {
            init_rayleigh_taylor(vbx, gbx,
                                 pres.array(mfi),
                                 vel.array(mfi),
                                 rho.array(mfi),
                                 trac.array(mfi),
                                 domain, dx, problo, probhi);
        }
        else if (12 == m_probtype)
        {
            init_periodic_tracer(vbx, gbx,
                                 pres.array(mfi),
                                 vel.array(mfi),
                                 rho.array(mfi),
                                 trac.array(mfi),
                                 domain, dx, problo, probhi);
        }
        else if (21 == m_probtype)
        {
            init_double_shear_layer(vbx, gbx,
                                    pres.array(mfi),
                                    vel.array(mfi),
                                    rho.array(mfi),
                                    trac.array(mfi),
                                    domain, dx, problo, probhi);
        }
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

void incflo::init_taylor_green (Box const& vbx, Box const& /* gbx */,
                                Array4<Real> const& /* p */,
                                Array4<Real> const& vel,
                                Array4<Real> const& /* density */,
                                Array4<Real> const& /* tracer */,
                                Box const& /* domain */,
                                GpuArray<Real, AMREX_SPACEDIM> const& dx,
                                GpuArray<Real, AMREX_SPACEDIM> const& /* problo */,
                                GpuArray<Real, AMREX_SPACEDIM> const& /* probhi */)
{
    amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        Real x = (i+0.5)*dx[0];
        Real y = (j+0.5)*dx[1];
        vel(i,j,k,0) =  std::sin(amr_wind::utils::two_pi()*x) * std::cos(amr_wind::utils::two_pi()*y);
        vel(i,j,k,1) = -std::cos(amr_wind::utils::two_pi()*x) * std::sin(amr_wind::utils::two_pi()*y);
        vel(i,j,k,2) = 0.0;
    });
}

void incflo::init_taylor_green3d (Box const& vbx, Box const& /* gbx */,
                                  Array4<Real> const& /* p */,
                                  Array4<Real> const& vel,
                                  Array4<Real> const& /* density */,
                                  Array4<Real> const& /* tracer */,
                                  Box const& /* domain */,
                                  GpuArray<Real, AMREX_SPACEDIM> const& dx,
                                  GpuArray<Real, AMREX_SPACEDIM> const& /* problo */,
                                  GpuArray<Real, AMREX_SPACEDIM> const& /* probhi */)
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

void incflo::init_couette (Box const& vbx, Box const& /* gbx */,
                           Array4<Real> const& /* p */,
                           Array4<Real> const& vel,
                           Array4<Real> const& /* density */,
                           Array4<Real> const& /* tracer  */,
                           Box const& domain,
                           GpuArray<Real, AMREX_SPACEDIM> const& /* dx */,
                           GpuArray<Real, AMREX_SPACEDIM> const& /* problo */,
                           GpuArray<Real, AMREX_SPACEDIM> const& /* probhi */)
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


void incflo::init_rayleigh_taylor (Box const& vbx, Box const& /* gbx */,
                                   Array4<Real> const& /* p */,
                                   Array4<Real> const& vel,
                                   Array4<Real> const& density,
                                   Array4<Real> const& /* tracer */,
                                   Box const& /* domain */,
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

void incflo::init_periodic_tracer (Box const& vbx, Box const& /* gbx */,
                                   Array4<Real> const& /* p */,
                                   Array4<Real> const& vel,
                                   Array4<Real> const& /* density */,
                                   Array4<Real> const& tracer,
                                   Box const& /* domain */,
                                   GpuArray<Real, AMREX_SPACEDIM> const& dx,
                                   GpuArray<Real, AMREX_SPACEDIM> const& problo,
                                   GpuArray<Real, AMREX_SPACEDIM> const& probhi)
{
    Real L = probhi[0]-problo[0];
    Real C = amr_wind::utils::two_pi() / L;
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

void incflo::init_double_shear_layer (Box const& vbx, Box const& /* gbx */,
                                      Array4<Real> const& /* p */,
                                      Array4<Real> const& vel,
                                      Array4<Real> const& /* density */,
                                      Array4<Real> const& tracer,
                                      Box const& /* domain */,
                                      GpuArray<Real, AMREX_SPACEDIM> const& dx,
                                      GpuArray<Real, AMREX_SPACEDIM> const& /* problo */,
                                      GpuArray<Real, AMREX_SPACEDIM> const& /* probhi */)
{

    amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        Real x = (i+0.5) * dx[0];
        Real y = (j+0.5) * dx[1];
        vel(i,j,k,0) = std::tanh(30.0*(0.25-std::abs(y-0.5)));
        vel(i,j,k,1) = 0.05*std::sin(amr_wind::utils::two_pi()*x);
        vel(i,j,k,2) = 0.0;

        Real r = std::sqrt((x-0.5)*(x-0.5) + (y-0.25)*(y-0.25));
        if (r < .1)
            tracer(i,j,k,0) = 0.0;
        else
            tracer(i,j,k,0) = 0.01;
    });


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
