#ifndef INCFLO_PROB_I_H_
#define INCFLO_PROB_I_H_

// This header is include by incflo.H.  These functions are members of incflo class.

    void init_plane_poiseuille (amrex::Box const& vbx, amrex::Box const& gbx,
                                amrex::Array4<amrex::Real> const& p,
                                amrex::Array4<amrex::Real> const& vel,
                                amrex::Array4<amrex::Real> const& density,
                                amrex::Array4<amrex::Real> const& tracer,
                                amrex::Box const& domain,
                                amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& dx,
                                amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& problo,
                                amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& probhi);

#endif
