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

    for (MFIter mfi(ld.density); mfi.isValid(); ++mfi)
    {
        const Box& vbx = mfi.validbox();
        const Box& gbx = mfi.fabbox();
        if (32 == probtype) 
        {
            init_plane_poiseuille(vbx, gbx,
                                  ld.p.array(mfi),
                                  ld.velocity.array(mfi),
                                  ld.density.array(mfi),
                                  ld.tracer.array(mfi),
                                  ld.eta.array(mfi),
                                  domain, dx, problo, probhi);
        }
        else
        {
            amrex::Abort("prob_init_fluid: unknown probtype");
        };
    }
}

void incflo::init_plane_poiseuille (amrex::Box const& vbx, amrex::Box const& gbx,
                                    amrex::Array4<amrex::Real> const& p,
                                    amrex::Array4<amrex::Real> const& vel,
                                    amrex::Array4<amrex::Real> const& density,
                                    amrex::Array4<amrex::Real> const& tracer,
                                    amrex::Array4<amrex::Real> const& eta,
                                    amrex::Box const& domain,
                                    GpuArray<Real, AMREX_SPACEDIM> const& dx,
                                    GpuArray<Real, AMREX_SPACEDIM> const& problo,
                                    GpuArray<Real, AMREX_SPACEDIM> const& probhi)
{
    Real dzinv = 1.0 / domain.length(2);
    const auto dlo = amrex::lbound(domain);
    const auto dhi = amrex::ubound(domain);
    Real lmu = this->mu;
    Real lrho = this->ro_0;
    if (32 == probtype)
    {
        Real v = this->ic_v;
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

            eta(i,j,k) = lmu;
        });
    }
    else
    {
        amrex::Abort("Unknown plane poiseuille probtype");
    };
}

