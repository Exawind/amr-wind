#ifndef PROB_BC_H_
#define PROB_BC_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_PhysBCFunct.H>

struct IncfloVelFill
{
    int probtype;
    amrex::GpuArray<amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>, AMREX_SPACEDIM*2> bcv;

    AMREX_GPU_HOST
    constexpr IncfloVelFill (int a_probtype,
                             amrex::GpuArray<amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>,
                                             AMREX_SPACEDIM*2> const& a_bcv)
        : probtype(a_probtype), bcv(a_bcv) {}

    AMREX_GPU_DEVICE
    void operator() (const amrex::IntVect& iv, amrex::Array4<amrex::Real> const& vel,
                     const int dcomp, const int numcomp,
                     amrex::GeometryData const& geom, const amrex::Real time,
                     const amrex::BCRec* bcr, const int bcomp,
                     const int orig_comp) const
    {
        using namespace amrex;

        // do something for external Dirichlet (BCType::ext_dir)
        const int i = iv[0];
        const int j = iv[1];
        const int k = iv[2];

        const Box& domain_box = geom.Domain();

        for (int vcomp = 0; vcomp < AMREX_SPACEDIM; ++vcomp)
        {
            const BCRec& bc = bcr[bcomp+vcomp];

            if (bc.lo(0) == BCType::ext_dir and i < domain_box.smallEnd(0))
            {
                vel(i,j,k,vcomp) = bcv[Orientation(Direction::x,Orientation::low)][vcomp];

                if (31 == probtype and vcomp == 0)
                {
                    Real y = (j+0.5)*(1.0/domain_box.length(1));
                    vel(i,j,k,0) *= 6. * y * (1.-y);
                } else
                if (311 == probtype and vcomp == 0)
                {
                    Real z = (k+0.5)*(1.0/domain_box.length(2));
                    vel(i,j,k,0) *= 6. * z * (1.-z);
                } else
                if (41 == probtype and vcomp == 0)
                {
                    Real z = (k+0.5)*(1.0/domain_box.length(2));
                    vel(i,j,k,0) =  0.5 * z;
                }
            }
            else if (bc.hi(0) == BCType::ext_dir and i > domain_box.bigEnd(0))
            {
                vel(i,j,k,vcomp) = bcv[Orientation(Direction::x,Orientation::high)][vcomp];
            }

            if (bc.lo(1) == BCType::ext_dir and j < domain_box.smallEnd(1))
            {
                vel(i,j,k,vcomp) = bcv[Orientation(Direction::y,Orientation::low)][vcomp];
                if (32 == probtype and vcomp == 1)
                {
                    Real z = (k+0.5)*(1.0/domain_box.length(2));
                    vel(i,j,k,1) *= 6. * z * (1.-z); 
                }
                if (322 == probtype and vcomp == 1)
                {
                    Real x = (i+0.5)*(1.0/domain_box.length(0));
                    vel(i,j,k,1) *= 6. * x * (1.-x); 
                }
            }
            else if (bc.hi(1) == BCType::ext_dir and j > domain_box.bigEnd(1))
            {
                vel(i,j,k,vcomp) = bcv[Orientation(Direction::y,Orientation::high)][vcomp];
            }

            if (bc.lo(2) == BCType::ext_dir and k < domain_box.smallEnd(2))
            {
                vel(i,j,k,vcomp) = bcv[Orientation(Direction::z,Orientation::low)][vcomp];
                if (33 == probtype and vcomp == 2)
                {
                    Real x = (i+0.5)*(1.0/domain_box.length(0));
                    vel(i,j,k,2) *= 6. * x * (1.-x); 
                }
                if (333 == probtype and vcomp == 2)
                {
                    Real y = (j+0.5)*(1.0/domain_box.length(1));
                    vel(i,j,k,2) *= 6. * y * (1.-y); 
                }
            }
            else if (bc.hi(2) == BCType::ext_dir and k > domain_box.bigEnd(2))
            {
                vel(i,j,k,vcomp) = bcv[Orientation(Direction::z,Orientation::high)][vcomp];
            }
        }
    }
};

struct IncfloDenFill
{
    int probtype;
    amrex::GpuArray<amrex::Real,AMREX_SPACEDIM*2> bcv;

    AMREX_GPU_HOST
    constexpr IncfloDenFill (int a_probtype, amrex::GpuArray<amrex::Real, AMREX_SPACEDIM*2> const& a_bcv)
        : probtype(a_probtype), bcv(a_bcv) {}

    AMREX_GPU_DEVICE
    void operator() (const amrex::IntVect& iv, amrex::Array4<amrex::Real> const& rho,
                     const int dcomp, const int numcomp,
                     amrex::GeometryData const& geom, const amrex::Real time,
                     const amrex::BCRec* bcr, const int bcomp,
                     const int orig_comp) const
    {
        using namespace amrex;

        // do something for external Dirichlet (BCType::ext_dir)
        const int i = iv[0];
        const int j = iv[1];
        const int k = iv[2];

        const Box& domain_box = geom.Domain();
        const BCRec& bc = bcr[bcomp];

        if (bc.lo(0) == BCType::ext_dir and i < domain_box.smallEnd(0))
        {
            rho(i,j,k) = bcv[Orientation(Direction::x,Orientation::low)];
        }
        else if (bc.hi(0) == BCType::ext_dir and i > domain_box.bigEnd(0))
        {
            rho(i,j,k) = bcv[Orientation(Direction::x,Orientation::high)];
        }

        if (bc.lo(1) == BCType::ext_dir and j < domain_box.smallEnd(1))
        {
            rho(i,j,k) = bcv[Orientation(Direction::y,Orientation::low)];
        }
        else if (bc.hi(1) == BCType::ext_dir and j > domain_box.bigEnd(1))
        {
            rho(i,j,k) = bcv[Orientation(Direction::y,Orientation::high)];
        }

        if (bc.lo(2) == BCType::ext_dir and k < domain_box.smallEnd(2))
        {
            rho(i,j,k) = bcv[Orientation(Direction::z,Orientation::low)];
        }
        else if (bc.hi(2) == BCType::ext_dir and k > domain_box.bigEnd(2))
        {
            rho(i,j,k) = bcv[Orientation(Direction::z,Orientation::high)];
        }
    }
};

struct IncfloTracFill
{
    int probtype, ntrac;
    amrex::GpuArray<amrex::Real const*, AMREX_SPACEDIM*2> bcv;

    AMREX_GPU_HOST
    constexpr IncfloTracFill (int a_probtype, int a_ntrac,
                              amrex::GpuArray<amrex::Real const*, AMREX_SPACEDIM*2> const& a_bcv)
        : probtype(a_probtype), ntrac(a_ntrac), bcv(a_bcv) {}

    AMREX_GPU_DEVICE
    void operator() (const amrex::IntVect& iv, amrex::Array4<amrex::Real> const& tracer,
                     const int dcomp, const int numcomp,
                     amrex::GeometryData const& geom, const amrex::Real time,
                     const amrex::BCRec* bcr, const int bcomp,
                     const int orig_comp) const
    {
        using namespace amrex;

        // do something for external Dirichlet (BCType::ext_dir)
        const int i = iv[0];
        const int j = iv[1];
        const int k = iv[2];

        const Box& domain_box = geom.Domain();

        for (int n = 0; n < ntrac; ++n)
        {
            const BCRec& bc = bcr[bcomp+n];

            if (bc.lo(0) == BCType::ext_dir and i < domain_box.smallEnd(0))
            {
                tracer(i,j,k,n) = bcv[Orientation(Direction::x,Orientation::low)][n];
            }
            else if (bc.hi(0) == BCType::ext_dir and i > domain_box.bigEnd(0))
            {
                tracer(i,j,k,n) = bcv[Orientation(Direction::x,Orientation::high)][n];
            }

            if (bc.lo(1) == BCType::ext_dir and j < domain_box.smallEnd(1))
            {
                tracer(i,j,k,n) = bcv[Orientation(Direction::y,Orientation::low)][n];
            }
            else if (bc.hi(1) == BCType::ext_dir and j > domain_box.bigEnd(1))
            {
                tracer(i,j,k,n) = bcv[Orientation(Direction::y,Orientation::high)][n];
            }

            if (bc.lo(2) == BCType::ext_dir and k < domain_box.smallEnd(2))
            {
                tracer(i,j,k,n) = bcv[Orientation(Direction::z,Orientation::low)][n];
            }
            else if (bc.hi(2) == BCType::ext_dir and k > domain_box.bigEnd(2))
            {
                tracer(i,j,k,n) = bcv[Orientation(Direction::z,Orientation::high)][n];
            }
        }
    }
};

struct IncfloForFill
{
    int probtype;

    AMREX_GPU_HOST
    constexpr IncfloForFill (int a_probtype) : probtype(a_probtype) {}

    AMREX_GPU_DEVICE
    void operator() (const amrex::IntVect& iv, amrex::Array4<amrex::Real> const& vel,
                     const int dcomp, const int numcomp,
                     amrex::GeometryData const& geom, const amrex::Real time,
                     const amrex::BCRec* bcr, const int bcomp,
                     const int orig_comp) const
    {}
};

#endif
